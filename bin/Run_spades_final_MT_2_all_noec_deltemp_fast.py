#!/usr/bin/env python
import pdb
#pdb.set_trace()
import glob
import os
from argparse import ArgumentParser
from multiprocessing import Pool,cpu_count,active_children,Manager
from subprocess import Popen
from Concatenate_contigs_all_v4_extend_for_HCbk import *
import multiprocessing
import pickle
import time
from datetime import timedelta
import pysam

def format_time(seconds):
    """Return human-readable time (min or hr)."""
    if seconds < 60:
        return f"{seconds :.1f} seconds"
    elif seconds < 3600:
        return f"{seconds / 60:.1f} min"
    else:
        return f"{seconds / 3600:.1f} hr"


script_path = os.path.dirname(os.path.abspath( __file__ ))
code_path = script_path + "/" 
parser = ArgumentParser(description="run local assembly by spades:")
# parser.add_argument('--bam_file','-bam',help="Required parameter; BAM file, called by bwa mem",required=True)
parser.add_argument('--sample_fastq','-fq', help="wgs fastq")
parser.add_argument('--num_threads','-nt', help="number of threads",default=30)
parser.add_argument('--num_threads_spades','-nt_spades', help="number of threads for spades",default=5)
parser.add_argument('--out_dir','-o', help="out dir")
parser.add_argument('--sample','-sp', help="sample")
parser.add_argument('--minicontig_dir','-minicontig_o', help="minicontig dir")
parser.add_argument('--chr_start','-start',type=int,help="chromosome start from", default=1)
parser.add_argument('--chr_end','-end',type=int,help="chromosome end by", default=23)
parser.add_argument('--read_type','-r',help = 'read type, single-end or paired-end. please input "SE" or "PE"', choices=['SE','PE'])
args = parser.parse_args()
num_threads_spades = int(args.num_threads_spades)

#new func
# def extract_read_by_idx_PE(f, out_file, read_idxs):
#     fw = open(out_file, 'w')
#     for idx in read_idxs[::2]:
#         f.seek(idx)
#         temp_read = list()
#         for _ in range(8):
#             line = f.readline()
#             fw.write(line)
#             temp_read.append(line)
#         r1_name = temp_read[0].rstrip("\n").split()[0]
#         r2_name = temp_read[4].rstrip("\n").split()[0]
#         assert r1_name == r2_name, "read names not match!"
#     fw.close()
#     return 

def reverse_complement(sequence):
    """
    Returns the reverse complement of a DNA sequence.
    
    Parameters:
        sequence (str): The input DNA sequence (uppercase or lowercase).
        
    Returns:
        str: The reverse complement of the DNA sequence.
    """
    complement = str.maketrans("ACGTacgt", "TGCAtgca")
    return sequence.translate(complement)[::-1]

def get_seq_from_bam(bam, idx, ):
    bam.seek(idx)
    read = next(bam)
    qname = read.qname
    # if qname!= exp_qname:
    #     print(f"{qname} and exp qname {exp_qname} don't match at {idx}")
    #     exit()
    if read.is_reverse:
        seq = reverse_complement(read.seq)
        qual = read.qual[::-1]
    else:
        seq = read.seq 
        qual = read.qual
    return qname,seq,qual

def extract_read_by_idx_PE(f,  out_file, read_idxs, ):
    fw = open(out_file, 'w')
    # filename = os.path.basename(out_file).split('_')[0]
    for pos in sorted(read_idxs)[::2]:
        idx = pos
        f.seek(idx)
        temp_read = list()
        for _ in range(8):
            line = f.readline()
            fw.write(line)
            temp_read.append(line)
        try:
            r1_name = temp_read[0].rstrip("\n").split()[0]
            r2_name = temp_read[4].rstrip("\n").split()[0]
        except:
            print('meow',idx,temp_read[0], temp_read[1])
            exit()
        # try:
        assert r1_name == r2_name, f"read names not match!{r1_name} {r2_name} {pos}"
        # except:
        #     print(r1_name, r2_name, pos)
        #     exit()
    fw.close()
    return 
    
def extract_read_by_idx_SE(f, out_file, read_idxs):
    fw = open(out_file, 'w')
    for idx in read_idxs:
        f.seek(idx)
        for _ in range(4):
            line = f.readline()
            fw.write(line)
    fw.close()
    return 

#new func
def extract_reads_from_sample_fastq(sample_fastq,one_file_fastq, read_idxs, read_type="PE"):
    # start_time = time.time()
    with open(sample_fastq,"r") as f:
        if read_type == "PE":
            extract_read_by_idx_PE(f,  one_file_fastq, read_idxs)
        else:
            extract_read_by_idx_SE(f, one_file_fastq, read_idxs)
    return




def use_spades(sample_fastq,one_file_fastq,read_idxs, fq_time_log,log_file, out_dir, read_type):
    start_time = time.time()
    extract_reads_from_sample_fastq(sample_fastq, one_file_fastq, read_idxs,  read_type)
    end_time1 = time.time()
    if read_type == 'PE':
        rt_flag = '--12'
    else:
        rt_flag = '-s'

    if os.path.exists(out_dir):
        os.system("rm -r " + out_dir)
    try:
        use_cmd = code_path + "SPAdes-3.13.0-Linux/bin/" + "spades.py -t " + str(num_threads_spades) + " --only-assembler %s "%rt_flag + one_file_fastq + " -o " + out_dir 
    except:
        use_cmd = "spades.py -t " + str(num_threads_spades) + " --only-assembler %s "%rt_flag + one_file_fastq + " -o " + out_dir 
    # with open('command.log','a') as f:
    #     f.write(use_cmd)
    Popen(use_cmd,shell=True).wait()
    end_time2 = time.time()
    # lock = threading.Lock()
    # with lock:
    #     with open(log_file,'a') as fw:
    #         fw.write(os.path.basename(one_file_fastq)+'\n')
    elapsed1 = format_time(end_time1 - start_time)
    elapsed2 = format_time(end_time2 - end_time1)
    elapsed3 = format_time(end_time2 - start_time)
    cmd = f"echo {elapsed1} {elapsed2} {elapsed3} >> {fq_time_log}"
    os.system(cmd)

    use_cmd = f"echo {os.path.basename(one_file_fastq)} >> {log_file}"
    print(use_cmd)
    Popen(use_cmd,shell=True).wait()



# def del_one_dir_1(out_dir):
#     temp_dir = out_dir + "_temp/"
#     use_cmd_1 = "mv " + out_dir + " " + temp_dir
#     Popen(use_cmd_1,shell=True).wait()
#     use_cmd_2 = "mkdir " + out_dir
#     Popen(use_cmd_2,shell=True).wait()
#     contig_file = temp_dir + "contigs.fasta"
#     if os.path.exists(contig_file):
#         use_cmd_3 = "mv " + contig_file + " " + out_dir
#         Popen(use_cmd_3,shell=True).wait()
#     rm_cmd = "rm -rf " + temp_dir
#     Popen(rm_cmd,shell=True).wait()


def del_dir(out_dir):
    rm_cmd = "rm -rf " + out_dir
    Popen(rm_cmd,shell=True).wait()

def del_file(file_dir):
    rm_cmd = "rm -f " + file_dir
    Popen(rm_cmd,shell=True).wait()

# TODO: here the read pickle file and sample fastq file were hard coded according to CAN's test files
# need to further modify them to be more general. ALSO, this is only for chr1 test, do not use this script for other chromosomes directly
def run_spades_all(sample_fastq, chr_start,chr_end,output_dir,num_of_threads,minicontig_dir,read_type,sample):
    import logging
    ## set logger
    logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')
    logger = logging.getLogger(" ")
    if not os.path.exists(minicontig_dir):
        os.system("mkdir -p " + minicontig_dir)


    # hard coded read pickle file and sample fastq file
    # sample_fastq = "/data/maiziezhou_lab/CanLuo/Aquila_stLFR2/Experiments/V350158928/chr1/chr1_fastq/HG002_chr1.fq"
    # read_idx_pickle = "/data/maiziezhou_lab/CanLuo/Aquila_stLFR2/Experiments/V350158928/chr1/Local_Assembly_by_chunks/chr1_files_cutPBHC/HP_Reads_merged.p"
    # chr_start=1
    # chr_end=1
    ##################################################
    #num_of_threads = multiprocessing.cpu_count() - 10
    # parent_dir = os.path.dirname(os.path.normpath(output_dir))
    error_file = f"{output_dir}/error_log"
    for chr_num in range(chr_start, chr_end+1):
        # sample_fastq = f"{parent_dir}/chr{chr_num}_fastq/{sample}_chr{chr_num}.fq"
        in_dir = output_dir + "chr" + str(chr_num) + "_files_cutPBHC/"
        log_file = in_dir+"/assmebly.log"
        if os.path.exists(log_file):
            os.system("rm -r " + log_file)
        else:
            os.system("touch " + log_file)

        read_idx_pickle = f"{in_dir}/HP_Reads_merged.p"
        output_file = "Aquila_cutPBHC_minicontig" + "_chr" + str(chr_num) + ".fasta" # this file contains all minicontigs for this chromosome
        minicontig_file = minicontig_dir+output_file
        if os.path.exists(minicontig_file):
            os.system("rm " + minicontig_file)
        ### hard coded part
        with open(read_idx_pickle, 'rb') as handle:
            read_idx_dict = pickle.load(handle)

        # total_num = len(read_idx_dict)
        
        spades_batch = list()
        async_results = None
        cmplt_dirs = list()
        cmplt_fastqs = list()

        n_total = len(read_idx_dict)
        time_log_file = in_dir + "/process.log"

        if os.path.exists(time_log_file):
            os.system("rm -r " + time_log_file)
        else:
            os.system("touch " + time_log_file)

        start_all = time.time()
        cmplt_cnt = 0
        prev_asm = glob.glob(in_dir +"/fastq_by_*_spades_assembly")
        print(prev_asm)
        [del_dir(asm) for asm in prev_asm]

        prev_asm = glob.glob(in_dir +"/fastq_by_*_hp*.fastq")
        print(prev_asm)
        [del_file(asm) for asm in prev_asm]
        # exit()
        contig_num = 1
        fq_time_log = in_dir +"/fq_time.log"
        with open(fq_time_log,'w') as f_fq_t:
            f_fq_t.write("extract_fq spades total\n")
        with Pool(num_of_threads) as pool:
            print(list(read_idx_dict.keys())[0][0])
            for key, read_idxs in read_idx_dict.items():
                s, e, hp = key
                one_file_fastq = in_dir + "/fastq_by_" + str(s) + "_" + str(e) + "_" + hp + ".fastq"
                # extract_reads_from_sample_fastq(bamfile,sample_fastq, one_file_fastq, read_idxs, error_file,fq_time_log, read_type)
                # with open(one_file_fastq, "w") as fw:
                #     fw.write(records)
        ### end of hard coded part
                out_dir = one_file_fastq[:-6] + "_spades_assembly"
                # spades_contig_file = out_dir + "/" + "contigs.fasta"

                with open(in_dir + "/cpu_count",'w') as f_cpu:
                    n = int(os.getenv("SLURM_CPUS_ON_NODE", 1))
                    
                    f_cpu.write(f"num of cpu: {os.cpu_count()}\n")
                    f_cpu.write(f"Allowed CPUs seen by Slurm: {n}\n")
                    
                # logger.info("==============start running")
                # use_spades(sample_fastq,one_file_fastq,read_idxs, fq_time_log,log_file, out_dir, read_type)
                # exit()

                if len(spades_batch) < num_of_threads:
                    # spades_batch.append((one_file_fastq, log_file, out_dir, read_type))
                    #=============
                    #     WARNING !!!
                    #   If you update below line, remember to update the cmplt_dirs and cmplt_fastqs
                    #    otherwise you will delete important files which you don't want!!!!
                    #=============
                    spades_batch.append((sample_fastq,one_file_fastq,read_idxs, fq_time_log,log_file, out_dir, read_type))
                else:
                    if async_results is not None: # wait for previous batch to finish
                        async_results.wait()
                        # collect results from completed dirs
                        contig_num = Concatenate_start(in_dir,minicontig_dir,output_file,contig_num) # note: move to append mode in Concatenate_start function
                        # remove completed dirs from disk
                        pool.starmap(del_dir, [(cmplt_dir,) for cmplt_dir in cmplt_dirs])
                        pool.starmap(del_file, [(one_file_fastq,) for one_file_fastq in cmplt_fastqs])
                        cmplt_cnt+=len(cmplt_fastqs)
                        cur_time = time.time()
                        total_elapsed = cur_time - start_all
                        avg_time = total_elapsed / cmplt_cnt
                        remaining = avg_time * (n_total - cmplt_cnt)

                        used_fmt = format_time(total_elapsed)
                        remain_fmt = format_time(remaining)

                        msg = (f"[{cmplt_cnt}/{n_total}] done in {used_fmt}, ETA {remain_fmt}\n")
                        print(msg, time_log_file)
                        with open(time_log_file, 'a') as fw:
                            fw.write(msg)

                        cmplt_dirs = list() # reset completed dirs list just in case
                        cmplt_fastqs = list() # reset completed fastq list just in case
                    
                    batch_input = spades_batch
                    cmplt_dirs = [x[-2] for x in spades_batch] # store batch output dirs
                    cmplt_fastqs = [x[1] for x in spades_batch] # store batch fastq files

                    async_results = pool.starmap_async(use_spades, batch_input) # submit current batch
                    # spades_batch = [(one_file_fastq, log_file, out_dir, read_type)] # start new batch
                    #=============
                    #     WARNING !!!
                    #   If you update below line, remember to update the cmplt_dirs and cmplt_fastqs
                    #    otherwise you will delete important files which you don't want!!!!
                    #=============
                    spades_batch = [(sample_fastq,one_file_fastq,read_idxs, fq_time_log,log_file, out_dir, read_type)]
            # process remaining jobs
            if spades_batch:
                if async_results is not None:
                    async_results.wait()
                    # collect results from completed dirs
                    contig_num = Concatenate_start(in_dir,minicontig_dir,output_file,contig_num) # note: move to append mode in Concatenate_start function
                    # remove completed dirs from disk
                    pool.starmap(del_dir, [(cmplt_dir,) for cmplt_dir in cmplt_dirs])
                    pool.starmap(del_file, [(one_file_fastq,) for one_file_fastq in cmplt_fastqs])
                    cmplt_dirs = list() # reset completed dirs list just in case
                    cmplt_fastqs = list() # reset completed fastq list just in case

                batch_input = spades_batch
                cmplt_dirs = [x[-2] for x in spades_batch] # store batch output dirs
                cmplt_fastqs = [x[1] for x in spades_batch] # store batch fastq files
                async_results = pool.starmap_async(use_spades, batch_input) # submit current batch
                async_results.wait()
                # final collection of results
                contig_num = Concatenate_start(in_dir,minicontig_dir,output_file,contig_num) # note: move to append mode in Concatenate_start function
                # remove completed dirs from disk
                pool.starmap(del_dir, [(cmplt_dir,) for cmplt_dir in cmplt_dirs])
                pool.starmap(del_file, [(one_file_fastq,) for one_file_fastq in cmplt_fastqs])
                cmplt_dirs = list() # reset completed dirs list just in case
                cmplt_fastqs = list() # reset completed fastq list just in case
       
    print("All Done~")

# backup
# def run_spades_all(chr_start,chr_end,output_dir,num_of_threads,minicontig_dir,read_type):
#     #num_of_threads = multiprocessing.cpu_count() - 10
#     for chr_num in range(chr_start, chr_end+1):
#         in_dir = output_dir + "chr" + str(chr_num) + "_files_cutPBHC/" 
#         count = 1
#         fastq_files_all = sorted(glob.glob(in_dir +  "fastq_by_*.fastq"),key=os.path.getsize,reverse=True)
#         total_num = len(fastq_files_all)
#         pool = Pool(num_of_threads)
#         out_dir_list = []
#         for one_file_fastq in fastq_files_all:
#             one_file = one_file_fastq[:-6]
#             out_dir = one_file + "_spades_assembly"
#             spades_contig_file = out_dir + "/" + "contigs.fasta"
#             if os.path.exists(spades_contig_file):
#                 count += 1
#                 #print("using existing " + spades_contig_file)
#                 out_dir_list.append(out_dir)
#             else:
#                 count += 1
#                 pool.apply_async(use_spades,(one_file_fastq,out_dir,read_type))
#                 out_dir_list.append(out_dir)
#             if (count - 1)%num_of_threads == 0 or (count - 1) == total_num:
#                 pool.close()
#                 while len(active_children()) > 1:
#                     time.sleep(0.5)
#                 pool.join()
                
#                 # delete temp files
#                 """
#                 if out_dir_list != []:
#                     pool = Pool(num_of_threads)
#                     for one_dir in out_dir_list:
#                         pool.apply_async(del_one_dir_1,(one_dir,"xin"))
#                     pool.close()
#                     while len(active_children()) > 1:
#                         time.sleep(0.5)
#                     pool.join()
#                 """

#                 if (count - 1) == total_num:
#                     print("finished chr" + str(chr_num))
#                 else:
#                     pool = Pool(num_of_threads)
#                     out_dir_list = []
#         output_file = "Aquila_cutPBHC_minicontig" + "_chr" + str(chr_num) + ".fasta"
#         Concatenate_start(in_dir,minicontig_dir,output_file,"xin")   
        
#         # delete assembly files
        
#         time.sleep(5)
#         pool = Pool(num_of_threads)
#         count = 1
#         for one_file_fastq in fastq_files_all:
#             one_file = one_file_fastq[:-6]
#             out_dir = one_file + "_spades_assembly"
#             count += 1
#             pool.apply_async(del_dir,(out_dir,"xin"))
#             if (count - 1)%num_of_threads == 0 or (count - 1) == total_num:
#                 pool.close()
#                 while len(active_children()) > 1:
#                     time.sleep(0.5)
#                 pool.join()
#                 if (count - 1) == total_num:
#                     print("finished chr" + str(chr_num))
#                 else:
#                     pool = Pool(num_of_threads)
#         time.sleep(5)
       
#     print("All Done~")






if __name__ == "__main__":
    # bamfile = args.bam_file
    sample_fastq = args.sample_fastq
    output_dir = args.out_dir
    minicontig_dir = args.minicontig_dir
    chr_start = args.chr_start
    chr_end = args.chr_end
    num_of_threads = int(args.num_threads)
    read_type = args.read_type
    sample = args.sample
    run_spades_all( sample_fastq, chr_start,chr_end,output_dir,num_of_threads,minicontig_dir,read_type,sample)

