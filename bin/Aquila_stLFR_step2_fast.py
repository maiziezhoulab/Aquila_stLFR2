#!/usr/bin/env python
import pdb
#pdb.set_trace()
from subprocess import Popen
from argparse import ArgumentParser
import os
import sys
import pickle
import os.path
script_path = os.path.dirname(os.path.abspath( __file__ ))
code_path = script_path + "/" 
__author__ = "Xin Zhou@Stanford"
parser = ArgumentParser(description="Author: xzhou15@cs.stanford.edu\n",usage='use "python3 %(prog)s --help" for more information')
# parser.add_argument('--bam_file','-bam',help="Required parameter; BAM file, called by bwa mem",required=True)
# parser.add_argument('--fastq_file','-fq', help="wgs fastq")
parser.add_argument('--chr_number','-chr',type=int,help="chromosome number, eg. 1,2,3...22")
# parser.add_argument('--chr_start','-start',type=int,help="chromosome start from, default = 1", default=1)
# parser.add_argument('--chr_end','-end',type=int,help="chromosome end by, default = 23", default=23)
parser.add_argument('--out_dir','-o', help="Directory to store assembly results, default = ./Assembly_results",default="./Asssembly_results")
parser.add_argument('--reference','-ref', help="Required parameter; reference fasta file",required=True)
parser.add_argument('--num_threads','-t',type=int,help="number of threads, default = 30, this correponds to number of small files get assembled simulateoulsy", default=30)
parser.add_argument('--num_threads_spades','-t_spades',type=int,help="number of threads for spades, default = 5", default=5)
parser.add_argument('--block_len_use','-bl',type=int,help="phase block len threshold, default = 100000",default=100000)
parser.add_argument('--read_type','-r',help = 'read type, single-end or paired-end. please input "SE" or "PE"', choices=['SE','PE'], default = 'PE')
parser.add_argument('--sample_name','-name',help="Required parameter; sample name you can define, for example, S12878",required=True)
args = parser.parse_args()

def read_ref(fasta_file,chr_num,out_dir):
    f = open(fasta_file,"r")
    count = 0
    ref_seq = ""
    for line in f:
        if count > 0:
            data = line.rsplit()
            ref_seq += data[0]
        count += 1
    print("total_len for chr" + str(chr_num))
    print(len(ref_seq))
    pickle.dump(ref_seq, open(out_dir + "ref_seq_chr" + str(chr_num) +  ".p","wb"))


def extract_ref_chr(ref_file,chr_num,out_dir):
    fw = open(out_dir + "genome_ref_chr" + str(chr_num) + ".fasta","w")
    f = open(ref_file,"r")
    flag = 0
    total_len = 0
    for line in f:
        data = line.rsplit()
        if data[0] == ">chr" + str(chr_num):
            fw.writelines(">" + str(chr_num) + "\n")
            flag = 1
        elif data[0][0] == ">" and flag == 1:
            break
        else:
            if flag == 1:
                total_len += len(data[0])
                fw.writelines(data[0] + "\n")
    print("chr" + str(chr_num) + ":")
    print(total_len)



def local_assembly_for_small_chunks( sample_fastq, chr_start,chr_end,num_threads,num_threads_spades,Local_Assembly_dir,Assembly_Contigs_dir,read_type, sample):
    use_cmd = "python " + code_path + "Run_spades_final_MT_2_all_noec_deltemp_fast.py" + " --num_threads " + str(num_threads) + " --num_threads_spades " + str(num_threads_spades) + " --chr_start " + str(chr_start) + " --chr_end " + str(chr_end) + " --out_dir " + Local_Assembly_dir + " --minicontig_dir " + Assembly_Contigs_dir +\
    " --read_type "+read_type + " --sample " + sample +  " --sample_fastq " + sample_fastq
    Popen(use_cmd,shell=True).wait()


def Complete_contiguity(chr_start,chr_end,Assembly_Contigs_dir,phase_blocks_cut_highconf_dir,cut_threshold,ref_file):
    use_cmd = "python " + code_path + "Concatenate_contigs_from_microcontigs_v2.py" + " --chr_start "  + str(chr_start) + " --chr_end " + str(chr_end) + " --out_dir " + Assembly_Contigs_dir + " --phase_cut_folder " + phase_blocks_cut_highconf_dir + " --cut_threshold " + str(cut_threshold) + " --ref_file " + ref_file 
    Popen(use_cmd,shell=True).wait()
  


def main():
    if len(sys.argv) == 1:
        Popen("python " + "Aquila_step2.py -h",shell=True).wait()
    else:
        # bam_file = args.bam_file
        # sample_fastq = args.fastq_file
        chr_start = args.chr_number
        chr_end = args.chr_number
        ref_file = args.reference
        cut_threshold = args.block_len_use
        num_threads = int(args.num_threads)
        num_threads_spades = int(args.num_threads_spades)
        read_type = args.read_type
        sample = args.sample_name

        
        Assembly_Contigs_dir = args.out_dir + "/Assembly_Contigs_files/"
        phase_blocks_cut_highconf_dir = args.out_dir + "/phase_blocks_cut_highconf/"
        Local_Assembly_dir = args.out_dir + "/Local_Assembly_by_chunks/"
        sample_fastq = f"{args.out_dir}/chr{chr_start}_fastq/{sample}_chr{chr_start}.fq"
        local_assembly_for_small_chunks(sample_fastq, chr_start,chr_end,num_threads,num_threads_spades,Local_Assembly_dir,Assembly_Contigs_dir,read_type,sample)
        Complete_contiguity(chr_start,chr_end,Assembly_Contigs_dir,phase_blocks_cut_highconf_dir,cut_threshold,ref_file)
    
if __name__ == "__main__":
    main()
