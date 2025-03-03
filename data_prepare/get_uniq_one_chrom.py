import argparse
from argparse import ArgumentParser
parser = ArgumentParser(description="",
	usage='use "python3 %(prog)s --help" for more information',
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--input_path','-i')
parser.add_argument('--output_path','-o')
#parser.add_argument('--n_thread','-t', type = int, default = 22 )
#parser.add_argument('--delete_temp_file','-d', action='store_true')
args = parser.parse_args()
input_path = args.input_path
output_path = args.output_path
#n_thread = args.n_thread


with open(input_path,'r') as f:
	s = list(set(f.readlines()))


with open(output_path,'w') as f:
	f.writelines(s)
