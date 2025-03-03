from argparse import ArgumentParser
from subprocess import Popen
import os
parser = ArgumentParser(description="",usage='use "python3 %(prog)s --help" for more information')
# parser.add_argument('--read_bam_file','-rbam')
parser.add_argument('--contig_path','-contig')
parser.add_argument('--reference_path','-ref')
# parser.add_argument('--signature_dir','-sigd')
parser.add_argument('--output_dir','-o')
parser.add_argument('--chr_num','-chr', type = int)
parser.add_argument('--n_thread','-t', type = int, default = 10)
parser.add_argument('--header_file','-header',default = "/data/maiziezhou_lab/CanLuo/long_reads_project/DipPAV_pipeline/side_bin/bin_v1/vcf_header",
	help = "default = /data/maiziezhou_lab/CanLuo/long_reads_project/DipPAV_pipeline/side_bin/bin_v1/vcf_header")

args = parser.parse_args()
# read_bam_file = args.read_bam_file
contig_path = args.contig_path
reference_path = args.reference_path
# signature_dir = args.signature_dir
header_file = args.header_file
output_dir = args.output_dir
chr_num = args.chr_num
n_thread = args.n_thread

import logging
## set logger
logging.basicConfig(
format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(" ")
import os
global code_dir
code_dir=os.path.dirname(os.path.realpath(__file__))+'/'

if not os.path.exists(output_dir):
	os.makedirs(output_dir)


logger.info("refromat fasta...")
def reformat_fasta(infa,outfa):
	with open(outfa,'w') as f:
		with open(infa,'r') as fin:
			cnt = 0
			for line in fin:
				if line[0]=='>':
					cnt+=1
					if 'hp1' in line:
						line = '>contig_%d_hp1\n'%cnt 
					else:
						line = '>contig_%d_hp2\n'%cnt 
				f.write(line)
	return 

contig_rf = output_dir+'/contigs.fa'
reformat_fasta(contig_path, contig_rf )

contig_path = contig_rf

logger.info("align contig to reference...")

cmd ="python3 "+code_dir+"/align_contig_to_reference.py \
-i  %s \
-o  %s -chr %d --ref %s -t %d"%(
	contig_path,
	output_dir,
	chr_num,
	reference_path,
	n_thread)
Popen(cmd,shell=True).wait()


logger.info("extract signatures from contigs bam file...")
prefix = contig_path.split('/')[-1].split('.')[0]
bam_path = output_dir+'/'+prefix+'.sorted.bam'
cmd = "python3 "+code_dir+"/extract_contig_signature.py \
-bam %s -chr %d -contig %s -header %s -ref %s -o %s"%(bam_path,chr_num,
 contig_path,
 header_file,
 reference_path,
 output_dir)

Popen(cmd,shell=True).wait()


# logger.info("extract signatures from reads bam file...")

# cmd = "python3 "+code_dir+"/extract_reads_signature.py \
# -i %s -chr %d  -o %s"%(read_bam_file ,chr_num,
#  output_dir)
# Popen(cmd,shell=True).wait()


# logger.info("Filter out false positive...")
vcf_path = output_dir+"/dippav_variant_chr%d.vcf"%chr_num
# signature_dir = output_dir+'/reads_signature/'
# vcf_path_filtered = output_dir+"/dippav_variant_chr%d_filtered.vcf"%chr_num
# cmd = "python3 "+code_dir+"/FP_filter_v1.py \
# -i %s -sigd %s -o %s"%(vcf_path,signature_dir, vcf_path_filtered)
# Popen(cmd,shell=True).wait()



logger.info("Remove redundancy...")
# vcf_path_noredun = output_dir+"/dippav_variant_chr%d_no_redundancy.vcf"%chr_num
final_vcf_dir = output_dir+'/final_vcf/'
cmd = "python3 "+code_dir+"/remove_redundancy.py -i %s -o %s "%(
	vcf_path, final_vcf_dir)
Popen(cmd,shell=True).wait()


logger.info("truvari evaluation...")
eval_dir = final_vcf_dir+'eval/'
cmd = code_dir + "/truvari2_eval.sh %d %s %s dippav_variant_no_redundancy 200 0.1 0.1 0 0 "%(
	chr_num, final_vcf_dir, eval_dir)
Popen(cmd,shell=True).wait()




