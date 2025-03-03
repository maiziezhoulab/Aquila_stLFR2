from argparse import ArgumentParser
from subprocess import Popen
import os
parser = ArgumentParser(description="",usage='use "python3 %(prog)s --help" for more information')
parser.add_argument('--input_path','-i')
parser.add_argument('--output_dir','-o')
parser.add_argument('--chr_num','-chr',type=int)
parser.add_argument('--ref')
parser.add_argument('--n_thread','-t',type = int, default = 20)
args = parser.parse_args()
input_path = args.input_path
output_dir = args.output_dir
chr_num = args.chr_num
ref_path = args.ref
n_thread = args.n_thread
if not os.path.exists(output_dir):
	os.makedirs(output_dir)

prefix = input_path.split('/')[-1].split('.')[0]


def fasta_to_bam(fasta_path,output_dir,ref_path):
	file_name = fasta_path.split('/')[-1].split('.')[0]
	cmd = "minimap2 -a -x asm5 --cs -r2k -t %d \
	   %s \
	   %s \
	   > %s/%s.sam"%(n_thread,ref_path,fasta_path,output_dir,file_name)
	Popen(cmd,shell=True).wait()

	cmd = "samtools sort -m4G -@4 \
	   -o %s/%s.sorted.bam \
	   %s/%s.sam"%(output_dir,file_name,output_dir,file_name)
	Popen(cmd,shell=True).wait()

	cmd = "samtools index %s/%s.sorted.bam"%(output_dir,file_name)
	Popen(cmd,shell=True).wait()

	os.system("rm %s/%s.sam"%(output_dir,file_name))

	return output_dir+'/'+file_name+'.sorted.bam'


bam_file = fasta_to_bam(input_path,output_dir,ref_path)

# ## split fasta file to 2 haplotypes

# with open(input_path,'r') as f:
# 	sf = f.readlines()

# hps = []
# idxs = []
# seqs = []
# for i in range(len(sf)):
# 	line = sf[i]
# 	if "hp" in line:
# 		hps.append(line)
# 		idxs.append(i)

# for j in range(len(idxs)-1):
# 	seqs.append(sf[idxs[j]+1:idxs[j+1]])
# seqs.append(sf[idxs[-1]+1:-1])


# hp1lines = []
# hp2lines = []

# for i in range(len(hps)):
# 	if "hp1" in hps[i]:
# 		hp1lines.append(hps[i])
# 		hp1lines.extend(seqs[i])
# 	else:
# 		hp2lines.append(hps[i])
# 		hp2lines.extend(seqs[i])

# with open("%s/%s_hp1.fa"%( output_dir, prefix),'w') as f:
# 	f.writelines(hp1lines)
# with open("%s/%s_hp2.fa"%( output_dir, prefix),'w') as f:
# 	f.writelines(hp2lines)


# ## align to reference

# fasta_path_hp1 = output_dir+'/'+prefix+'_hp1.fa'
# fasta_path_hp2 = output_dir+'/'+prefix+'_hp2.fa'
# bam_hp1 = fasta_to_bam(fasta_path_hp1,output_dir,ref_path)
# bam_hp2 = fasta_to_bam(fasta_path_hp2,output_dir,ref_path)







