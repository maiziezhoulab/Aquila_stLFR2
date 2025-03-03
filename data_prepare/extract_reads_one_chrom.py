
import argparse
from argparse import ArgumentParser
parser = ArgumentParser(description="",
	usage='use "python3 %(prog)s --help" for more information',
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--fq_file','-fq')
parser.add_argument('--qn_file','-qn')
parser.add_argument('--output_dir','-o')

parser.add_argument('--out_prefix','-px', default = "NA24385")
parser.add_argument('--n_thread','-t', type = int, default = 5 )
parser.add_argument('--n_chunk','-nck', type = int, default = 100 )
parser.add_argument('--delete_temp_file','-d', action='store_true')
args = parser.parse_args()
fq_file = args.fq_file
qn_file = args.qn_file
out_prefix = args.out_prefix
output_dir = args.output_dir
n_thread = args.n_thread
n_chunk = args.n_chunk


import os
from tqdm import tqdm
import pandas as pd

from joblib import Parallel, delayed
from tqdm import tqdm
from subprocess import Popen


def find_real_bnd(file,fake_bt):
	''' This will return  
	the starting byte of the closest read name line 
	that is after you given byte

	Note:
	if your byte is already in the last line of the file, or exceeding the file size
	this function will return the ENDing byte of the file 
	'''
	filesize = os.path.getsize(file)
	if fake_bt >= filesize:
		print(f'''Given byte {fake_bt} is greater than or equal to the file size {filesize}
		return the file size {filesize} as the real cutting byte''' )
		return filesize

	f = open(file,'r')
	f.seek(fake_bt)
	while True:
		real_bt = f.tell()
		l = f.readline()
		if l:
			if (l[0]=='@') & ('BX:Z:' in l):
				break
		else:
			break
			
	f.close() 
	return real_bt

def cut_file_into_chunks(file, n_ck, output_dir):
	filesize = os.path.getsize(file)
	size_per_ck = int(filesize / n_ck) + 1
	assert size_per_ck * n_ck >= filesize 
	print(f"num chunks: {n_ck}")
	print(f"chunk size: {size_per_ck / (1024)**3} GB")

	fk_rgl = [  (i*size_per_ck, (i+1)* size_per_ck) for i in range(n_ck-1)]
	st = (n_ck-1) * size_per_ck 
	ed =  filesize 
	if st < ed:
		fk_rgl.append((st,ed))
	# print(fk_rgl)
	real_bnd_list = [  find_real_bnd(file,rg[1])   for rg in fk_rgl[:-1]]

	real_rgl = [ [0, real_bnd_list[0]]]

	for i in range(1,len(real_bnd_list)):
		real_rgl.append( [real_bnd_list[i-1], real_bnd_list[i]])
	
	real_rgl.append([real_bnd_list[-1], filesize])
	assert len(real_rgl) == len(fk_rgl)
	# print(real_rgl)

	df = pd.DataFrame(real_rgl, columns= ['start','end'])
	df.to_csv(output_dir+'/chunks.csv', index = False)
	print(df)



	return real_rgl







def process_one_part(fq_file, qn_file, out_dir, rg, prefix):
	st,ed = rg
	outfile  = out_dir+f'/{prefix}_{st}_{ed}.fq'
	print(f"loading qn file {qn_file}...")
	with open(qn_file,'r') as f:
		qnames = set(f.read().split('\n')[:-1])
	print(f"finished loading qn file {qn_file}")
	i=0
	with open(outfile,'w') as fw:
		f = open(fq_file,'r') 

		f.seek(st)
		while line := f.readline():
			i+=1
			if i%4==1:
				in_chrom_flag = 0
				qname = line[1:].split()[0].split('#')[0]
				if qname in qnames:
					in_chrom_flag = 1
					qnames.remove(qname)
			if in_chrom_flag:
				fw.write(line)
			cur_bt = f.tell()
			if cur_bt >= ed:
				break 
	del qnames
	return 




def reduce_files(outfiles,  final_file):

	cmd = f" cat {outfiles[0]} >  {final_file}"
	print(cmd)
	Popen(cmd, shell = True).wait()


	for file in tqdm(outfiles[1:], desc = "reduce files"):
		cmd = "cat " + file + " >> " + final_file
		Popen(cmd, shell = True).wait()

		if args.delete_temp_file:
			cmd = "rm " + file
			# print(cmd)
			Popen(cmd, shell = True).wait()
	return 



if not os.path.exists(output_dir):
	os.system("mkdir -p " + output_dir)

#------------- cut file into chunks

real_rgl = cut_file_into_chunks(fq_file, n_chunk, output_dir)


#------------- process by chunk
# sequences = Parallel(n_jobs=n_thread)(delayed(process_one_part)\
# 									  (fq_file, qn_file, output_dir, rg, out_prefix) for rg in tqdm(real_rgl))

outfiles = [output_dir+f'/{out_prefix}_{rg[0]}_{rg[1]}.fq' for rg in real_rgl]
final_file = output_dir + f"/merged_{out_prefix}.fq"
reduce_files(outfiles,  final_file)


