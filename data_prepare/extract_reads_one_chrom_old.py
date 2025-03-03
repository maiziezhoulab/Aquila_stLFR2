
import argparse
from argparse import ArgumentParser
parser = ArgumentParser(description="",
	usage='use "python3 %(prog)s --help" for more information',
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--fq_file_pattern','-fq')
parser.add_argument('--output_dir','-o')
parser.add_argument('--out_prefix','-px', default = "NA24385")
parser.add_argument('--qn_file','-qn')
parser.add_argument('--n_thread','-t', type = int, default = 22 )
parser.add_argument('--delete_temp_file','-d', action='store_true')
args = parser.parse_args()
fq_file_pattern = args.fq_file_pattern
qn_file = args.qn_file
out_prefix = args.out_prefix
output_dir = args.output_dir
n_thread = args.n_thread




def process_one_part(fq_file, qn_file, outfile):


	with open(qn_file,'r') as f:
		qnames = set(f.read().split('\n')[:-1])

	i=0

	with open(outfile,'w') as fw:
		with open(fq_file,'r') as f:
			for line in f:
				i+=1
				if i%8==1:
					in_chrom_flag = 0
					qname = line[1:].split()[0]
					if qname in qnames:
						in_chrom_flag = 1
						qnames.remove(qname)

				if in_chrom_flag:
					fw.write(line)

				# if (i!=1) and (i%10000000 == 0):
				# 	print(f"processed {i} lines")

	del qnames
	
	return 


from subprocess import Popen
def reduce_files(fq_file_pattern, output_dir, prefix):
	out_pattern = output_dir + '/'+fq_file_pattern.split('/')[-1]
	final_file = output_dir + '/' + prefix + '.fastq'

	cmd = "cat " + out_pattern + " > " + final_file
	print(cmd)
	Popen(cmd, shell = True).wait()

	if args.delete_temp_file:
		cmd = "rm " + out_pattern
		print(cmd)
		Popen(cmd, shell = True).wait()



	return 

import glob

fq_file_list = glob.glob(fq_file_pattern)
# print(fq_file_list)

import os
if not os.path.exists(output_dir):
	os.system("mkdir -p " + output_dir)

outfile_list = []
for fq_file in fq_file_list:
	prefix = fq_file.split('/')[-1]
	batch_file = output_dir+'/'+prefix
	outfile_list.append(batch_file)

# print(outfile_list)

from joblib import Parallel, delayed
from tqdm import tqdm
# n_thread = 10
sequences = Parallel(n_jobs=n_thread)(delayed(process_one_part)(fq_file_list[i], qn_file, outfile_list[i]) for i in tqdm(range(len(fq_file_list))))




reduce_files(fq_file_pattern, output_dir, out_prefix)

