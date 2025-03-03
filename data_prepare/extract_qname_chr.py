import argparse
from argparse import ArgumentParser
parser = ArgumentParser(description="",
	usage='use "python3 %(prog)s --help" for more information',
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--bam_file','-bam')
parser.add_argument('--out_dir','-o')
parser.add_argument('--chr_num','-chr', type = int, help = "1-22; optional, if not provided, run chr1-22")

parser.add_argument('--bin_size','-bin', type = int, default = 10000000)
parser.add_argument('--n_thread','-t', type = int, default = 22 )

args = parser.parse_args()
bam_file = args.bam_file
out_dir = args.out_dir
chr_num = args.chr_num
bin_size = args.bin_size
n_thread = args.n_thread




import pysam

import os

def cut_reference_into_bins(reference_lengths, bin_size):
    reference_bins = []

    for reference_name, reference_length in reference_lengths.items():
        bins = []
        for start in range(0, reference_length, bin_size):
            end = min(start + bin_size, reference_length)
            bins.append((f'{reference_name}:{start+1}-{end}'))
        reference_bins.extend(bins)
    print('number of bins:',len(reference_bins))
    print(reference_bins)

    return reference_bins


def get_ref_len(bamfile):
	
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


def extract_readnames_to_file(bam_file, region, outdir):
    # Ensure the output directory exists
    os.makedirs(outdir, exist_ok=True)
    
    # Initialize the BAM file reader
    
    bam = pysam.AlignmentFile(bam_file, "rb")

    # Fetch reads for the specified region
    try:
    	reads = bam.fetch(region=region)
    except:
    	print(region)
    	exit()

    # Create the output file path
    output_file = os.path.join(outdir, f"{region}.txt")

    # Open the output file for writing
    with open(output_file, "w") as output:
        # Iterate over reads and write read names to the output file
        for read in reads:
            output.write(f"{read.query_name}\n")

    # Close the BAM file
    bam.close()

    print(f"Read names extracted to {output_file}")


def merge_chunks(out_dir, chrom):

	
	cmd = f'cat {out_dir}/{chrom}:*.txt > {out_dir}/merged_{chrom}.txt'
	os.system(cmd)
	print(f"Qnames merged to {out_dir}/merged_{chrom}.txt")


if chr_num is None:
	chr_list = [ 'chr'+str(i+1) for i in range(22)]
else:
	assert 1<= chr_num <= 22
	chr_list = ['chr'+str(chr_num)]
	
reference_lengths = get_ref_len(bam_file)
reference_bins = cut_reference_into_bins(reference_lengths, bin_size)

#	ref_bin = [ b  for b in reference_bins if int(b.split(':')[0][3:]) == chr_num]
#	reference_bins = ref_bin 
print(chr_list)
print(len(reference_bins))
print(reference_bins)




from joblib import Parallel, delayed
from tqdm import tqdm

sequences = Parallel(n_jobs=n_thread)(delayed(extract_readnames_to_file)(bam_file, region, out_dir) for region in tqdm(reference_bins))



sequences = Parallel(n_jobs=n_thread)(delayed(merge_chunks)(out_dir, chrom) for chrom in tqdm(chr_list, desc = 'merge_chunks'))
