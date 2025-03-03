import argparse
from argparse import ArgumentParser
parser = ArgumentParser(description="",
	usage='use "python3 %(prog)s --help" for more information',
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--bam_file','-bam')
parser.add_argument('--out_dir','-o')
parser.add_argument('--chr_num','-chr', type = int, help = "1-22")
parser.add_argument('--fq_file_r1_L01','-fq1_L01')
parser.add_argument('--fq_file_r2_L01','-fq2_L01')
parser.add_argument('--fq_file_r1_L02','-fq1_L02')
parser.add_argument('--fq_file_r2_L02','-fq2_L02')
parser.add_argument('--fq_file_r1_L03','-fq1_L03')
parser.add_argument('--fq_file_r2_L03','-fq2_L03')
parser.add_argument('--fq_file_r1_L04','-fq1_L04')
parser.add_argument('--fq_file_r2_L04','-fq2_L04')
parser.add_argument('--index_dir_r1_L01','-index_r1_L01')
parser.add_argument('--index_dir_r2_L01','-index_r2_L01')
parser.add_argument('--index_dir_r1_L02','-index_r1_L02')
parser.add_argument('--index_dir_r2_L02','-index_r2_L02')
parser.add_argument('--index_dir_r1_L03','-index_r1_L03')
parser.add_argument('--index_dir_r2_L03','-index_r2_L03')
parser.add_argument('--index_dir_r1_L04','-index_r1_L04')
parser.add_argument('--index_dir_r2_L04','-index_r2_L04')


parser.add_argument('--bin_size','-bin', type = int, default = 10000000)
parser.add_argument('--n_thread','-t', type = int, default = 50 )
# parser.add_argument('--n_chunk','-nck', type = int, help = "n chunks for extracting reads from wgs fastq", default= 70)
parser.add_argument('--out_prefix','-px', default = "Sample")
args = parser.parse_args()
bam_file = args.bam_file
out_dir = args.out_dir
chr_num = args.chr_num
bin_size = args.bin_size
n_thread = args.n_thread
index_dir_r1_L01 = args.index_dir_r1_L01
index_dir_r2_L01 = args.index_dir_r2_L01
index_dir_r1_L02 = args.index_dir_r1_L02
index_dir_r2_L02 = args.index_dir_r2_L02
index_dir_r1_L03 = args.index_dir_r1_L03
index_dir_r2_L03 = args.index_dir_r2_L03
index_dir_r1_L04 = args.index_dir_r1_L04
index_dir_r2_L04 = args.index_dir_r2_L04
fqfile_r1_L01 = args.fq_file_r1_L01
fqfile_r2_L01 = args.fq_file_r2_L01
fqfile_r1_L02 = args.fq_file_r1_L02
fqfile_r2_L02 = args.fq_file_r2_L02
fqfile_r1_L03 = args.fq_file_r1_L03
fqfile_r2_L03 = args.fq_file_r2_L03
fqfile_r1_L04 = args.fq_file_r1_L04
fqfile_r2_L04 = args.fq_file_r2_L04
# n_ck = args.n_chunk
out_prefix = args.out_prefix



import pysam
from joblib import Parallel, delayed
from tqdm import tqdm
from subprocess import Popen
import os
code_dir = os.path.dirname(os.path.realpath(__file__))+'/'
import sys

# Add the directory path where your other script is located
sys.path.append(code_dir)


from retrieve_reads_multi_threads import retrieve_reads_multi_threads



def cut_reference_into_bins(reference_lengths, bin_size,out_dir):
    reference_bins = []

    for reference_name, reference_length in reference_lengths.items():
        bins = []
        for start in range(0, reference_length, bin_size):
            end = min(start + bin_size, reference_length)
            bins.append((f'{reference_name}:{start+1}-{end}'))
        reference_bins.extend(bins)
    print('number of bins:',len(reference_bins))
    print(reference_bins)

    outfile = out_dir+'/chunks.txt'
    with open(outfile,'w') as f:
        f.write('\n'.join(reference_bins)+'\n')

    return reference_bins


def get_ref_len(bamfile,chr_list):

    # Replace 'your_file.bam' with the path to your BAM or CRAM file
    bam_file = pysam.AlignmentFile(bamfile, 'rb')

    # Get reference sequence lengths
    reference_lengths = dict(zip(bam_file.references, bam_file.lengths))

    reference_lengths = {key:val for key,val in reference_lengths.items() if key in chr_list}

    # Print reference lengths
    for reference_name, reference_length in reference_lengths.items():
        print(f"Reference: {reference_name}, Length: {reference_length}")

    # Close the BAM or CRAM file when done
    bam_file.close()

    return reference_lengths


def extract_readnames_to_file(bam_file, region, outdir):
    # Ensure the output directory exists
    os.makedirs(outdir, exist_ok=True)

    # Initialize the BAM file reader
    print(region)
    bam = pysam.AlignmentFile(bam_file, "rb")

    # check readlen
    N = -1
    for read in bam.fetch(region = region):
        print(read.qname)
        cigar = read.cigar 
        if (cigar[0][0]!=5) and (cigar[-1][0]!=5):
            N = len(read.seq)
            break 
    print("read len:",N, region)
                 


    # Fetch reads for the specified region
    try:
        reads = bam.fetch(region=region)
    except:
        print(region)
        exit()

    # Create the output file path
    output_file = os.path.join(outdir, f"{region}.txt")
    used_file = os.path.join(outdir, f"{region}_used.txt")

    fq_file = os.path.join(outdir, f"{region}.fq")

    

    dc_r1 = {}
    dc_r2 = {}


    fw_used = open(used_file,'w')
    # Open the output file for writing
    with open(output_file, "w") as output:
        with open(fq_file, "w") as fw:
            # Iterate over reads and write read names to the output file
            for read in reads:
                cigar = read.cigar 
                if len(cigar) < 1:
                    # print(cigar)
                    continue
                
                if (cigar[0][0]!=5) and (cigar[-1][0]!=5):
                    if read.is_reverse:
                        seq = read.get_forward_sequence()
                        qual = read.qual[::-1]
                    else:
                        seq = read.seq 
                        qual = read.qual 
                    assert len(seq)==N
                    assert len(qual) ==N 
                    if read.is_read1:
                        dc_r1[read.qname] = [seq,qual]
                    elif read.is_read2:
                        dc_r2[read.qname] = [seq,qual]
                    qn = read.qname

                    if (read.qname in dc_r1) and (read.qname in dc_r2):
                        fw.write("@"+read.qname+'\n')
                        fw.write(dc_r1[qn][0]+'\n')
                        fw.write('+\n')
                        fw.write(dc_r1[qn][1]+'\n')
                        fw.write("@"+read.qname+'\n')
                        fw.write(dc_r2[qn][0]+'\n')
                        fw.write('+\n')
                        fw.write(dc_r2[qn][1]+'\n')
                        dc_r1.pop(qn)
                        dc_r2.pop(qn)
                        fw_used.write(f"{read.query_name}\n")
                                
                              
                output.write(f"{read.query_name}\n")

    # Close the BAM file
    bam.close()
    fw_used.close()

    print(f"Read names extracted to {output_file}")


def merge(in_list, out_file):
    cmd = 'cat '+ ' '.join(in_list) + " > " + out_file
    print(cmd)
    Popen(cmd, shell= True).wait()
    return 

def get_diff( all_name,used_name, outfile):
    with open(all_name,'r') as f:
        all = set(f.readlines())
    with open(used_name,'r') as f:
        used = set(f.readlines())
    unused = list(all - used)
    del all 
    del used 

    with open(outfile,'w') as f:
        f.writelines(unused)
    return 



def extract_fq_from_bam(chr_num, bam_file, out_dir, bin_size, n_thread):
    if chr_num is None:
        chr_list = [ 'chr'+str(i+1) for i in range(22)]
    else:
        assert 1<= chr_num <= 22
        chr_list = ['chr'+str(chr_num)]
        
    reference_lengths = get_ref_len(bam_file, chr_list)
    reference_bins = cut_reference_into_bins(reference_lengths, bin_size , out_dir)


    sequences = Parallel(n_jobs=n_thread)(delayed(extract_readnames_to_file)(bam_file, region, out_dir) for region in tqdm(reference_bins))

    all_fq = [out_dir+'/'+bin+'.fq' for bin in reference_bins]
    merged_fq = out_dir+'/merged.fq'
    merge(all_fq, merged_fq)


    all_name_list = [ out_dir+'/'+bin+'.txt' for bin in reference_bins]
    all_name  = out_dir+'/all_name.txt'
    merge(all_name_list, all_name)


    used_name_list = [ out_dir+'/'+bin+'_used.txt' for bin in reference_bins]
    used_name  = out_dir+'/used_name.txt'
    merge(used_name_list, used_name)


    unused_name_uniq  = out_dir+'/unused_name_uniq.txt'
    get_diff(all_name, used_name, unused_name_uniq)
    return 


def load_names(nf):
    with open(nf,'r') as f:
        names  = f.read().split('\n')[:-1]
    return names

def cut_chunks(names,n_ck, outdir, prefix):
    size_per_ck = int(len(names) / n_ck) + 1
    chunks  = [(i*size_per_ck, (i+1)*size_per_ck) for i in range(n_ck-1)]
    chunks.append((size_per_ck*(n_ck-1), len(names)))

    outfile = outdir+'/chunks_'+prefix+'.txt'
    with open(outfile,'w') as f:
        print(chunks, file = f)
    return chunks 

def merge(in_list, out_file):
    cmd = 'cat '+ ' '.join(in_list) + " > " + out_file
    print(cmd)
    Popen(cmd, shell= True).wait()
    return 

def split_qname_by_lane(qname_file):
    outfile_list =  [qname_file.replace(".txt","_L0%d.txt"%i) for i in range(1,5)]
    fw_list = []
    for ofile in outfile_list:
        fw = open(ofile, 'w')
        fw_list.append(fw)

    with open(qname_file,'r') as f:
        for line in f:
            qname = line[:-1]
            ln = int(qname[11])
            fw = fw_list[ln -1 ]
            fw.write(line)

    for fw in fw_list:
        fw.close()

    



def extract_fq_from_fq(qname_file, out_dir, fastq_file, index_dir, prefix,n_thread, ):
    # qname_file = out_dir+'/unused_name_uniq.txt'
    output_path = out_dir+f'/merged_{prefix}.fq'
    retrieve_reads_multi_threads(fastq_file, qname_file, output_path, index_dir, n_thread)

    return 

def merge_fastq(input1, input2, output):
    with open(input1, 'r') as file1, open(input2, 'r') as file2, open(output, 'w') as output_file:
        while True:
            try:
                for _ in range(4):
                    if _ == 0:
                        # line = 
                        output_file.write(next(file1).split()[0].split('#')[0]+'\n')
                    else:

                        output_file.write(next(file1))
                for _ in range(4):
                    if _ == 0:
                         output_file.write(next(file2).split()[0].split('#')[0]+'\n')
                    else:
                        output_file.write(next(file2))
            except StopIteration:
                break


def final_merge(out_dir):
    cmd = f"cat {out_dir}/merged.fq {out_dir}/merged_r12.fq > {out_dir}/merged_all.fq"
    print(cmd)
    Popen(cmd,shell = True,).wait()
    return 


def load_read_names_from_fastq(fastq_file):
    read_names = []
    i = 0
    with open(fastq_file, 'r') as f:
        for line in f:
            i+=1
            if i%4 == 1:
                # Extract the read name from the FASTQ header line
                read_name = line.strip().split()[0][1:]
                read_names.append(read_name)
                if i%8 ==  5 :
                    assert read_name ==  read_names[-1]
                    
    return set(read_names)

def load_read_names_from_file(read_name_file):
    read_names = []
    with open(read_name_file, 'r') as f:
        for line in f:
            read_names.append(line.strip())
    return set(read_names)

def sainity_check(out_dir):
    fastq_file = out_dir+'/merged_all.fq'
    read_name_file = out_dir+'/all_name.txt'
    uniq_read_names_rn_file = out_dir+'/uniq_read_names_rn.txt'
    uniq_read_names_fq_file = out_dir+'/uniq_read_names_fq.txt'

    fastq_read_names = load_read_names_from_fastq(fastq_file)
    read_name_file_read_names = load_read_names_from_file(read_name_file)

    print("uniq names in fq:",len(fastq_read_names))
    print("uniq names in rn:",len(read_name_file_read_names))

    uniq_read_names_rn = list(read_name_file_read_names - fastq_read_names)
    uniq_read_names_fq = list(fastq_read_names - read_name_file_read_names )

    print("in rn but not in fq:", len(uniq_read_names_rn))
    print("in fq but not in rn:", len(uniq_read_names_fq))

    with open(uniq_read_names_fq_file,'w') as f:
        s = '\n'.join(uniq_read_names_fq)+'\n'
        f.write(s)

    with open(uniq_read_names_rn_file,'w') as f:
        s = '\n'.join(uniq_read_names_rn)+'\n'
        f.write(s)


    assert len(uniq_read_names_fq) == 0



if not os.path.exists(out_dir):
    os.system("mkdir  -p " + out_dir)

# ##---------------- extract reads from bam file first (exclude non-paired, hard clipped reads)
# extract_fq_from_bam(chr_num, bam_file, out_dir, bin_size, n_thread, )


# ##------------------ extract reads from wgs fastq file using index (non-paired, hard-clipped in bam file)
qname_file = out_dir+'/unused_name_uniq.txt'
outfile_list =  [qname_file.replace(".txt","_L0%d.txt"%i) for i in range(1,5)]
split_qname_by_lane(qname_file)

extract_fq_from_fq(outfile_list[0], out_dir, fqfile_r1_L01, index_dir_r1_L01, 'r1_L01', n_thread)
extract_fq_from_fq(outfile_list[0], out_dir, fqfile_r2_L01, index_dir_r2_L01, 'r2_L01', n_thread)

extract_fq_from_fq(outfile_list[1], out_dir, fqfile_r1_L02, index_dir_r1_L02, 'r1_L02', n_thread)
extract_fq_from_fq(outfile_list[1], out_dir, fqfile_r2_L02, index_dir_r2_L02, 'r2_L02', n_thread)

extract_fq_from_fq(outfile_list[2], out_dir, fqfile_r1_L03, index_dir_r1_L03, 'r1_L03', n_thread)
extract_fq_from_fq(outfile_list[2], out_dir, fqfile_r2_L03, index_dir_r2_L03, 'r2_L03', n_thread)

extract_fq_from_fq(outfile_list[3], out_dir, fqfile_r1_L04, index_dir_r1_L04, 'r1_L04', n_thread)
extract_fq_from_fq(outfile_list[3], out_dir, fqfile_r2_L04, index_dir_r2_L04, 'r2_L04', n_thread)
    
r1_list =  [out_dir + "/merged_r1_L0%d.fq"%i for i in range(1,5)]
r1_fq = out_dir+"/merged_r1.fq"
cmd = "cat "+ " ".join(r1_list) + " > " + r1_fq
Popen(cmd, shell = True).wait()

r2_list =  [out_dir + "/merged_r2_L0%d.fq"%i for i in range(1,5)]
r2_fq = out_dir+"/merged_r2.fq"
cmd = "cat "+ " ".join(r2_list) + " > " + r2_fq
Popen(cmd, shell = True).wait()

# ##--------------------- pair reads file
merge_fastq(out_dir+'/merged_r1.fq', out_dir+'/merged_r2.fq', out_dir+'/merged_r12.fq')


# #------------ merge reads from different source
final_merge(out_dir)

# ------------ check if any reads are missing; if reads are interleaf
sainity_check(out_dir)