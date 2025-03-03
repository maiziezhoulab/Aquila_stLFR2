from argparse import ArgumentParser
from collections import defaultdict
from tqdm import tqdm
parser = ArgumentParser(description="",usage='use "python3 %(prog)s --help" for more information')
parser.add_argument('--input_path','-i')
parser.add_argument('--output_path','-o')
parser.add_argument('--reference','-ref')
parser.add_argument('--header_file','-header',default = "/data/maiziezhou_lab/CanLuo/Software/Aquila_stLFR/bin/vcf_header",
	help = "default = /data/maiziezhou_lab/CanLuo/Software/Aquila_stLFR/bin/vcf_header")

parser.add_argument('--delete_temp_file','-d', action='store_true')
args = parser.parse_args()
input_path = args.input_path
output_path = args.output_path
reference = args.reference
header_file = args.header_file


with open(input_path,'r') as f:
	lines = []
	for line in f:
		if ("SVTYPE=INS" in line) or ("SVTYPE=DEL" in line):
			lines.append(line)

with open(header_file,'r') as f:
	header = f.read()

with open(reference,'r') as f:
	dc = defaultdict(list)
	for line in f:
		if line[0]=='>':
			chrom = line[1:-1]
		else:
			dc[chrom].append(line[:-1])
	for chrom in dc:
		dc[chrom] = ''.join(dc[chrom])

with open(output_path,'w') as f:
	f.write(header)
	for line in tqdm(lines):
		data = line.split()
		chrom = data[0]
		pos = int(data[1])
		seq = dc[chrom]
		if "SVTYPE=DEL" in line:
			ref = seq[pos-1: pos + len(data[3])].upper()
			alt = seq[pos-1].upper()
			assert ref[1:].lower() == data[3].lower()
			data[3] = ref 
			data[4] = alt 
			data[7] = "SVTYPE=DEL;SVLEN=%d"%(len(alt)-len(ref))
		else:
			ref = seq[pos-1].upper()
			alt = (ref+data[4]).upper()
			data[3] = ref 
			data[4] = alt 
			data[7] = "SVTYPE=INS;SVLEN=%d"%(len(alt)-len(ref))
		line = '\t'.join(data)+'\n'
		f.write(line)








