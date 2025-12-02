#!/usr/bin/env python
import pdb
#pdb.set_trace()
from subprocess import Popen
from argparse import ArgumentParser
import os
import sys
script_path = os.path.dirname(os.path.abspath( __file__ ))
code_path = script_path + "/" 
__author__ = "Xin Zhou@Stanford"
parser = ArgumentParser(description="Author: can.luo@vanderbilt.edu\n",usage='use "python3 %(prog)s --help" for more information')
parser.add_argument('--fastq_file','-f',help="Required parameter; stLFR FASTQ file with paired reads",required=True)
parser.add_argument('--bam_file','-bam',help="Required parameter; BAM file, called by bwa mem",required=True)
parser.add_argument('--vcf_file','-v',help="Required parameter; VCF file, called by FreeBayes",required=True)
parser.add_argument('--chr_number','-chr',type=int,help="chromosome number, eg. 1,2,3...22")
parser.add_argument('--max_mem','-mem', help = "GB", default= 100, type = int)
# parser.add_argument('--chr_start','-start',type=int,help="chromosome start from, default = 1", default=1)
# parser.add_argument('--chr_end','-end',type=int,help="chromosome end by,default = 23", default=22)
parser.add_argument('--sample_name','-name',help="Required parameter; sample name you can define, for example, S12878",required=True)
parser.add_argument('--out_dir','-o', help="Directory to store assembly results, default = ./Assembly_results",default="./Asssembly_results")
parser.add_argument('--uniq_map_dir','-uniq_dir', help="Required Parameter; Directory for 100-mer uniqness, run ./install to download it",required=True)
parser.add_argument('--num_threads','-t_chr',type=int,help="number of threads, default = 8 (recommended)", default=8)
parser.add_argument('--block_threshold','-bt',type=int,help="phase block threshold, default = 200000",default=200000)
parser.add_argument('--block_len_use','-bl',type=int,help="phase block len threshold, default = 100000",default=100000)
parser.add_argument('--read_type','-r',help = 'read type, single-end or paired-end. please input "SE" or "PE"', choices=['SE','PE'])
args = parser.parse_args()

import logging
## set logger
logging.basicConfig(
format='%(asctime)s %(levelname)-8s %(message)s',
level=logging.INFO,
datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(" ")

from Generate_highconf_cut_profile import gen_highconf_profile

def Get_fragment_files(bam_file,vcf_file,chr_start,chr_end,h5_dir,num_threads,sample_name,read_type):
    use_cmd = "python3 " + code_path + "Run_h5_all_multithreads_v4.py" + " --bam_file " + bam_file + \
    " --vcf_file " + vcf_file + " --sample_name " + \
        sample_name + " --chr_num " \
            + str(chr_start) \
                + " --mbq 13 --mmq 20 --boundary 50000 " + \
                    " --num_threads " + str(num_threads) + " --out_dir " + h5_dir +\
                    " --read_type "+read_type
    print(use_cmd)
    logger.info("******************************************************\n\n")
    logger.info("               Generate molecule H5 file              \n\n")
    logger.info("******************************************************\n\n")
    Popen(use_cmd,shell=True).wait()


  
def Get_highconf_profile(bam_file,chr_start,chr_end,HighConf_file_dir,uniq_map_dir, threads):
    logger.info("******************************************************\n\n")
    logger.info("     Generate high confidence cut profile             \n\n")
    logger.info("******************************************************\n\n")
    gen_highconf_profile(bam_file,chr_start,chr_end,13,20,HighConf_file_dir,uniq_map_dir,threads)

def Haplotying_fragments(chr_start,chr_end,phased_file_dir,h5_dir,sample_name, threads):
    use_cmd = "python3 " + code_path + "Run_phase_alg_multithreads2.py" + " --chr_start " + str(chr_start) + " --chr_end " + str(chr_end) + \
        " --overlap_threshold 3 --support_threshold 5 " + " --out_dir "\
              + phased_file_dir  + " --h5_dir " + h5_dir + " --sample_name " + sample_name\
              + " -t " + str(threads)
    print(use_cmd)
    logger.info("******************************************************\n\n")
    logger.info("                 Phasing molecules                    \n\n")
    logger.info("******************************************************\n\n")
    Popen(use_cmd,shell=True).wait()


def Cut_phase_blocks(chr_start,chr_end,block_threshold,block_len_use,phase_blocks_cut_highconf_dir,phased_file_dir,HighConf_file_dir):
    use_cmd = "python3 " + code_path + "Cut_phaseblock_for_phased_h5_v4.0_highconf_v2.py" + " --chr_start " + str(chr_start) + " --chr_end " + str(chr_end) + " --block_threshold " + str(block_threshold) + " --block_len_use " + str(block_len_use)  + " --out_dir " +  phase_blocks_cut_highconf_dir + " --in_dir " + phased_file_dir + " --highconf_profile_dir " + HighConf_file_dir
    print(use_cmd)
    logger.info("******************************************************\n\n")
    logger.info("                   Cut phase blocks                   \n\n")
    logger.info("******************************************************\n\n")
    Popen(use_cmd,shell=True).wait()


#def Get_fastq_files_total(bam_file,chr_start,chr_end,num_threads,Raw_fastqs_dir,Sorted_bam_dir):
def Get_fastq_files_total(fastq_file,chr_start,chr_end,num_threads,Raw_fastqs_dir,h5_dir,sample_name):
    #use_cmd = "python3 " + code_path + "Read_fastqs_from_sortedbam_v2.py " + " --num_threads " + str(num_threads) + " --chr_start " + str(chr_start) + " --chr_end " + str(chr_end) + " --out_dir " + Raw_fastqs_dir + " --bam_file " + bam_file + " --bam_dir " + Sorted_bam_dir
    use_cmd = "python3 " + code_path + "Run_all_chunks.py " + " --out_dir_prefix " + " Raw_fastqs "  + " --h5_dir " + h5_dir + " --sample_name " + sample_name + " --num_threads " + str(num_threads) + " --chr_start "  + str(chr_start)  + " --chr_end " + str(chr_end) + " --out_dir " + Raw_fastqs_dir + " --fastq_file " + fastq_file
    print(use_cmd)
    logger.info("******************************************************\n\n")
    logger.info("               Split fastq by chromosome              \n\n")
    logger.info("******************************************************\n\n")
    Popen(use_cmd,shell=True).wait()


def Extract_reads_for_small_chunks_old(chr_start,chr_end,h5_dir,phase_blocks_cut_highconf_dir,Local_Assembly_dir,Raw_fastqs_dir,block_len_use,sample_name,num_threads,read_type):
    use_cmd = "python3 " + code_path + "Run_extract_reads_for_smallchunks_all_lessmem.py" + \
        " --phase_cut_folder " + phase_blocks_cut_highconf_dir + " --chr_start " + str(chr_start) +\
              " --chr_end " + str(chr_end) + " --out_dir " + Local_Assembly_dir + " --h5_folder " + h5_dir + \
                  " --fastq_folder " + Raw_fastqs_dir + " --cut_threshold " + str(block_len_use) + \
        " --sample_name " + sample_name + " --num_threads " +   str(num_threads) +\
        " --read_type "+read_type
    print(use_cmd)
    logger.info("******************************************************\n\n")
    logger.info("               Split fastq by phase block             \n\n")
    logger.info("******************************************************\n\n")
    Popen(use_cmd,shell=True).wait()



def Extract_reads_for_one_chromosome(bam_file, input_dir,fastq_file,chr_num, num_threads, sample_name):
    out_dir = f"{input_dir}/chr{chr_num}_fastq/"
    use_cmd = "python3 " + code_path + "BAM2FASTQ_By_Chromosome.py" + \
        " --bam_file " + bam_file + \
        " --fastq_file " + fastq_file + \
        " --prefix " + sample_name + " --chromosome chr" + str(chr_num) + \
        " -t " +  str(num_threads) +\
        " --out_dir "+ out_dir
    print(use_cmd)
    logger.info("******************************************************\n\n")
    logger.info("               Extract by chrom reads byte positions             \n\n")
    logger.info("******************************************************\n\n")
    Popen(use_cmd,shell=True).wait()

def Extract_reads_for_small_chunks(input_dir,chr_num, num_threads,max_mem, sample_name):
    fastq_file_one_chr = f"{input_dir}/chr{chr_num}_fastq/{sample_name}_chr{chr_num}.fq"
    use_cmd = "python3 " + code_path + "Extract_qname_from_phased_molecule_cut_phase_blocks_v6.py" + \
        " --indir " + input_dir + " --chr_fastq " + fastq_file_one_chr + \
        " --sample_name " + sample_name + " --chr_num " + str(chr_num) + " --n_thread " +  str(num_threads) +\
        " --max_mem "+str(max_mem)
    print(use_cmd)
    logger.info("******************************************************\n\n")
    logger.info("               Split fastq by phase block             \n\n")
    logger.info("******************************************************\n\n")
    Popen(use_cmd,shell=True).wait()

def clean_folder(out_dir, chr_num):
    chrom = 'chr'+str(chr_num)
    use_cmd = f"cd {out_dir}; rm -r {chrom}_fastq H5_for_molecules  HighConf_file results_phased_probmodel; " +\
            f"cd phase_blocks_cut_highconf; rm chr*_HC_breakpoint.bed chr*.phased_final_cut_by_100000  chr*.phased_final_cut_by_100000_phase_blocks* "
    logger.info("******************************************************\n\n")
    logger.info("               Clean up Step1 intermediate files             \n\n")
    logger.info("******************************************************\n\n")
    logger.info(use_cmd)
    Popen(use_cmd,shell=True).wait()

def main():
    if len(sys.argv) == 1:
        Popen("python3 " + "Aquila_step1.py -h",shell=True).wait()
    else:
        fastq_file = args.fastq_file    # add
        bam_file = args.bam_file
        vcf_file = args.vcf_file
        
        chr_start = args.chr_number
        chr_end = args.chr_number
        max_mem = args.max_mem
        block_len_use = args.block_len_use
        block_threshold = args.block_threshold
        uniq_map_dir = args.uniq_map_dir + "/"
        num_threads = int(args.num_threads)
        #num_threads_for_bwa_mem = int(args.num_threads_for_bwa_mem)
        sample_name = args.sample_name
        h5_dir = args.out_dir + "/H5_for_molecules/"
        HighConf_file_dir = args.out_dir + "/HighConf_file/"
        phased_file_dir = args.out_dir + "/results_phased_probmodel/"
        phase_blocks_cut_highconf_dir = args.out_dir + "/phase_blocks_cut_highconf/"
        Raw_fastqs_dir = args.out_dir + "/Raw_fastqs_chr" + str(chr_start) + "_" + str(chr_end) +  "/"
        Sorted_bam_dir = args.out_dir + "/sorted_bam/"
        Local_Assembly_dir = args.out_dir + "/Local_Assembly_by_chunks/"
        read_type = args.read_type

        assert read_type in ['SE','PE']


        Get_fragment_files(bam_file,vcf_file,chr_start,chr_end,h5_dir,num_threads,sample_name,read_type) 
        Get_highconf_profile(bam_file,chr_start,chr_end,HighConf_file_dir,uniq_map_dir, num_threads)
        Haplotying_fragments(chr_start,chr_end,phased_file_dir,h5_dir,sample_name, num_threads)
        Cut_phase_blocks(chr_start,chr_end,block_threshold,block_len_use,phase_blocks_cut_highconf_dir,phased_file_dir,HighConf_file_dir) 
        Extract_reads_for_one_chromosome(bam_file, args.out_dir,fastq_file,args.chr_number, num_threads, sample_name)
        Extract_reads_for_small_chunks(args.out_dir, args.chr_number, num_threads,max_mem, sample_name)
        # clean_folder(args.out_dir, args.chr_number)
        
if __name__ == "__main__":
    main()
