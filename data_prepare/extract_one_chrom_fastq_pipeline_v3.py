import argparse
from argparse import ArgumentParser
parser = ArgumentParser(description="",
	usage='use "python3 %(prog)s --help" for more information',
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--bam_file','-bam')
parser.add_argument('--out_fq_file','-o')
parser.add_argument('--temp_dir','-temp')
parser.add_argument('--chr_num','-chr', type = int, help = "1-22")
parser.add_argument('--qname_separator','-qsep', default= None, help = "qname field separator;optional, if not provided, will use anything before the first tab as qname in the qname line")
parser.add_argument('--lane_pos','-lp', type = int, help = "lane character position in read name; 0 based; optional, only needed if the fastqs are multiple lane")
parser.add_argument('--fq_file_r1','-fq1', nargs = "+")
parser.add_argument('--fq_file_r2','-fq2', nargs = "+")
parser.add_argument('--index_dir_r1','-index_r1', nargs = "+")
parser.add_argument('--index_dir_r2','-index_r2', nargs = "+")
# parser.add_argument('--bin_size','-bin', type = int, default = 10000000)
parser.add_argument('--n_thread','-t', type = int, default = 50 )
parser.add_argument('--out_prefix','-px', default = "Sample")
parser.add_argument('--debug','-d', action='store_true', help = "if set, will not delete intermediate files")
parser.add_argument('--strict','-s', action='store_true', help = "strict mode; if set, will perform sanity check in the end, this might take from 10mins to 30 mins depending on how large your file is")
args = parser.parse_args()
bam_file = args.bam_file
out_fq_file = args.out_fq_file
out_dir = args.temp_dir
chr_num = args.chr_num
# bin_size = args.bin_size
n_thread = args.n_thread
lane_pos = args.lane_pos
qname_separator = args.qname_separator
index_dir_r1_list = args.index_dir_r1
index_dir_r2_list = args.index_dir_r2


fqfile_r1_list = args.fq_file_r1
fqfile_r2_list = args.fq_file_r2


out_prefix = args.out_prefix
global debug 
debug = args.debug
strict = args.strict

print("qsep:", qname_separator)

import pysam
from joblib import Parallel, delayed
from tqdm import tqdm
from subprocess import Popen
import os
code_dir = os.path.dirname(os.path.realpath(__file__))+'/'
import sys

# Add the directory path where your other script is located
# sys.path.append(code_dir)


from retrieve_reads_multi_threads import retrieve_reads_multi_threads



def cut_reference_into_bins(reference_lengths, n_ck, out_dir,):
    reference_bins = []

    

    for reference_name, reference_length in reference_lengths.items():
        bins = []
        bin_size = int(reference_length/n_ck) + 1
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


def load_uniq(infile):
    with open(infile,'r') as f:
        s = set(f.readlines())
    return s 

def sanity_check_one_region(all_file, used_file, unused_file):
    s_all = load_uniq(all_file)
    s_used = load_uniq(used_file)
    s_unused = load_uniq(unused_file)
    s_diff = s_all -  (s_used|s_unused)  
    if len(s_diff):
        print( f" {used_file} and {unused_file} can not fully cover {all_file}, num uncovered is {len(s_diff)}")
        exit()
    
    assert len(s_diff)==0



def extract_readnames_to_file(bam_file, region, outdir):
    # Ensure the output directory exists
    os.makedirs(outdir, exist_ok=True)

    # Initialize the BAM file reader
    print(region)
    bam = pysam.AlignmentFile(bam_file, "rb")

    # check readlen
    N = -1
    for read in bam.fetch(region = region):
        # print(read.qname)
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
    unused_file = os.path.join(outdir, f"{region}_unused.txt")

    fq_file = os.path.join(outdir, f"{region}.fq")

    

    dc_r1 = {}
    dc_r2 = {}

    dc_check = {}

    fw_used = open(used_file,'w')
    # Open the output file for writing
    with open(output_file, "w") as output:
        with open(fq_file, "w") as fw:
            # Iterate over reads and write read names to the output file
            for read in reads:
                 

                output.write(f"{read.query_name}\n")
                dc_check[read.qname] = 1

                cigar = read.cigar
                if (len(cigar) >0) :                
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
                            # if used, remove qname from dc_check
                            dc_check.pop(qn)
                            fw_used.write(f"{read.query_name}\n")
                                
    unused_qname_list = list(dc_check.keys())


                              
            
    # Close the BAM file
    bam.close()
    fw_used.close()

    with open(unused_file,'w') as fw:
        if len(unused_qname_list):
            s = '\n'.join(unused_qname_list)+'\n'
        else:
            s = ''
        fw.write(s)

    sanity_check_one_region(output_file, used_file, unused_file)

    print(f"Read names extracted to {output_file}")




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

def write_uniq(infile, ofile):
    with open(infile, 'r') as f:
        s= list(set(f.readlines()))

    with open(ofile, 'w') as f:
        f.writelines(s)


def get_bins(chr_num, n_thread, bam_file):
    if chr_num is None:
        chr_list = [ 'chr'+str(i+1) for i in range(22)]
    else:
        assert 1<= chr_num <= 22
        chr_list = ['chr'+str(chr_num)]
        
    reference_lengths = get_ref_len(bam_file, chr_list)
    n_ck = n_thread
    reference_bins = cut_reference_into_bins(reference_lengths, n_ck , out_dir)
    return reference_bins

def extract_fq_from_bam(reference_bins,bam_file, out_dir,  n_thread):



    sequences = Parallel(n_jobs=n_thread)(delayed(extract_readnames_to_file)(bam_file, region, out_dir) for region in tqdm(reference_bins))

    # all_fq = [out_dir+'/'+bin+'.fq' for bin in reference_bins]
    # merged_fq = out_dir+'/merged.fq'
    # merge(all_fq, merged_fq)


    # all_name_list = [ out_dir+'/'+bin+'.txt' for bin in reference_bins]
    # all_name  = out_dir+'/all_name.txt'
    # merge(all_name_list, all_name)


    # used_name_list = [ out_dir+'/'+bin+'_used.txt' for bin in reference_bins]
    # used_name  = out_dir+'/used_name.txt'
    # merge(used_name_list, used_name)

    unused_name_list = [ out_dir+'/'+bin+'_unused.txt' for bin in reference_bins]
    unused_name  = out_dir+'/unused_name.txt'
    merge(unused_name_list, unused_name)


    unused_name_uniq  = out_dir+'/unused_name_uniq.txt'
    # get_diff(all_name, used_name, unused_name_uniq)
    write_uniq(unused_name, unused_name_uniq)
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

    if not debug:
        cmd = "rm " + ' '.join(in_list)
        Popen(cmd, shell= True).wait()
    return 

def split_qname_by_lane(qname_file, num_lane, lane_pos):
    outfile_list =  [qname_file.replace(".txt","_L0%d.txt"%i) for i in range(1,num_lane+1)]
    fw_list = []
    for ofile in outfile_list:
        fw = open(ofile, 'w')
        fw_list.append(fw)

    with open(qname_file,'r') as f:
        for line in f:
            qname = line[:-1]
            ln = int(qname[lane_pos])
            fw = fw_list[ln -1 ]
            fw.write(line)

    for fw in fw_list:
        fw.close()

    



def extract_fq_from_fq(qname_file, out_dir, fastq_file, index_dir, prefix,n_thread, qname_separator ):
    # qname_file = out_dir+'/unused_name_uniq.txt'
    output_path = out_dir+f'/merged_{prefix}.fq'
    retrieve_reads_multi_threads(fastq_file, qname_file, output_path, index_dir, n_thread, qname_separator)

    return 

def merge_fastq(input1, input2, output, qname_separator):
    with open(input1, 'r') as file1, open(input2, 'r') as file2, open(output, 'w') as output_file:
        while True:
            try:
                for _ in range(4):
                    if _ == 0:
                        # line = 
                        output_file.write(next(file1).split()[0].split(qname_separator)[0]+'\n')
                    else:

                        output_file.write(next(file1))
                for _ in range(4):
                    if _ == 0:
                         output_file.write(next(file2).split()[0].split(qname_separator)[0]+'\n')
                    else:
                        output_file.write(next(file2))
            except StopIteration:
                break


def final_merge(out_dir, all_fq_from_bam, out_fq_file):

    remove_dup(all_fq_from_bam + [f"{out_dir}/merged_r12.fq"], out_fq_file)
    fqs = ' '.join(all_fq_from_bam)

    cmd = f"cat {out_dir}/merged_r12.fq {fqs} > " + out_fq_file
    print(cmd)
    Popen(cmd,shell = True,).wait()

    if not debug:
        cmd = "rm "+ fqs 
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


def remove_dup(inlist, ofile):
    name_dc ={}
    dup_cnt = 0
    uniq_cnt = 0
    dup_log = ofile + ".dedup_log"
    with open(ofile,'w') as fw:
        for infile in inlist:
            with open(infile,'r') as f:
                i=0
                write_flag = 0
                for line in f:
                    i+=1
                    if i%8 ==1 :
                        qname = line[1:-1]
                        if qname not in name_dc:
                            name_dc[qname] =1
                            uniq_cnt+=1
                            write_flag = 1
                        else:
                            dup_cnt +=1
                    if write_flag:
                        fw.write(line)
                        
                    if i%8 == 0:
                        write_flag = 0
                    
    assert uniq_cnt == len(name_dc)
    print("DUP event: "+str(dup_cnt))
    print("Uniq event: "+str(uniq_cnt))
    with open(dup_log,'w') as fw:
        fw.write("DUP event: "+str(dup_cnt)+'\n')
        fw.write("Uniq event: "+str(uniq_cnt)+'\n')
    name_set = set(list(name_dc.keys()))
    return name_set


def extract_fq_from_fq_list(out_dir, fqfile_r1_list, fqfile_r2_list, index_dir_r1_list, index_dir_r2_list, n_thread, lane_pos, qname_separator):
    qname_file = out_dir+'/unused_name_uniq.txt'
    if len(fqfile_r1_list)>1:
        print("=============Multi lane mode")

        # check if input list is in order

        find_var_char(fqfile_r1_list)
        find_var_char(fqfile_r2_list)
        find_var_char(index_dir_r1_list)
        find_var_char(index_dir_r2_list)
        outfile_list =  [qname_file.replace(".txt","_L0%d.txt"%i) for i in range(1,len(fqfile_r1_list)+1)]
        # split_qname_by_lane(qname_file)
        num_lane = len(fqfile_r1_list)
        split_qname_by_lane(qname_file, num_lane, lane_pos)

        for i in range(len(fqfile_r1_list)):
            print("Lane %d ========="%i)
            print(outfile_list[i])
            print(fqfile_r1_list[i])
            print(fqfile_r2_list[i])
            print(index_dir_r1_list[i])
            print(index_dir_r2_list[i])


        for i in range(len(fqfile_r1_list)):
            extract_fq_from_fq(outfile_list[i], out_dir, fqfile_r1_list[i], index_dir_r1_list[i], 'r1_L0%d'%(i+1), n_thread, qname_separator)
            extract_fq_from_fq(outfile_list[i], out_dir, fqfile_r2_list[i], index_dir_r2_list[i], 'r2_L0%d'%(i+1), n_thread, qname_separator)
        
        r1_list =  [out_dir + "/merged_r1_L0%d.fq"%i for i in range(1,len(fqfile_r1_list)+1)]
        r1_fq = out_dir+"/merged_r1.fq"
        cmd = "cat "+ " ".join(r1_list) + " > " + r1_fq
        Popen(cmd, shell = True).wait()

        r2_list =  [out_dir + "/merged_r2_L0%d.fq"%i for i in range(1,len(fqfile_r1_list)+1)]
        r2_fq = out_dir+"/merged_r2.fq"
        cmd = "cat "+ " ".join(r2_list) + " > " + r2_fq
        Popen(cmd, shell = True).wait()

    else:
        print("=============Single lane mode")
        extract_fq_from_fq(qname_file, out_dir, fqfile_r1_list[0], index_dir_r1_list[0], 'r1', n_thread)
        extract_fq_from_fq(qname_file, out_dir, fqfile_r2_list[0], index_dir_r2_list[0], 'r2', n_thread)
    
    # ##--------------------- pair reads file
    merge_fastq(out_dir+'/merged_r1.fq', out_dir+'/merged_r2.fq', out_dir+'/merged_r12.fq', qname_separator)



def sanity_check(out_dir, fastq_file, reference_bins, fastq_read_names ):
    print("start final sanity check")

    all_name_list = [ out_dir+'/'+bin+'.txt' for bin in reference_bins]
    # all_name  = out_dir+'/all_name.txt'
    # merge(all_name_list, all_name)


    # fastq_file = out_dir+f'/chr{chr_num}.fq'
    read_name_file = out_dir+'/all_name.txt'
    if len(fastq_read_names)==0:
        fastq_read_names = load_read_names_from_fastq(fastq_file)

    read_name_file_read_names = []
    for read_name_file in tqdm(all_name_list, desc = "load all read names"):
        read_name_file_read_names.extend(list(load_read_names_from_file(read_name_file) ))
    read_name_file_read_names = set(read_name_file_read_names)

    print("uniq names in fq:",len(fastq_read_names))
    print("uniq names in rn:",len(read_name_file_read_names))

    diff_read_names = read_name_file_read_names - fastq_read_names
    print("in rn but not in fq:", len(diff_read_names))
    assert len(diff_read_names) == 0

    


def find_var_char(inlist):
    var_dc = {}
    ref  = [ str(i+1) for i in range(len(inlist))]
    for i in range(len(inlist[0])):

        chars = []
        for j in range(len(inlist)):
            chars.append(inlist[j][i])
        if len(set(chars))>1:
            var_dc[i] = chars 

    print(var_dc)
    
    var_idx = list(var_dc.keys())[0]

    assert len(var_dc)==1
    assert var_dc[var_idx] == ref
    print(var_dc)



if not os.path.exists(out_dir):
    os.system("mkdir  -p " + out_dir)




# ##---------------- extract reads from bam file first (exclude non-paired, hard clipped reads)
reference_bins = get_bins(chr_num, n_thread, bam_file)
all_fq_from_bam = [out_dir+'/'+bin+'.fq' for bin in reference_bins]
extract_fq_from_bam(reference_bins,  bam_file, out_dir,  n_thread, )


# ##------------------ extract reads from wgs fastq file using index (non-paired, hard-clipped in bam file)

extract_fq_from_fq_list(out_dir, fqfile_r1_list, fqfile_r2_list, 
                        index_dir_r1_list, index_dir_r2_list, 
                        n_thread, lane_pos,qname_separator)



# #------------ merge reads from different source
name_set_fq = {}
inlist = all_fq_from_bam + [f"{out_dir}/merged_r12.fq"]
name_set_fq = remove_dup(inlist , out_fq_file)

if strict:
    # ------------ check if any reads are missing; if reads are interleaf
    sanity_check(out_dir,out_fq_file, reference_bins,name_set_fq)

if not debug:
    #--------pass sanity check, delte all temp
    cmd = "rm -r " + out_dir
    Popen(cmd, shell = True).wait()