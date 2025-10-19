import pysam
from joblib import Parallel, delayed
from tqdm import tqdm
from subprocess import Popen
import os
import gzip
from multiprocessing import Pool
import tempfile
import argparse

def get_chromosome_length(bam_file, chromosome):
    """
    Get the length of a specific chromosome from a BAM file.
    Args:
        bam_file (str): Path to the BAM file.
        chromosome (str): Name of the chromosome (e.g., 'chr1').
    Returns:
        int: Length of the chromosome, or None if not found.
    """
    # Open the BAM file
    bam = pysam.AlignmentFile(bam_file, "rb")
    
    # Fetch the reference lengths and names
    references = bam.references
    lengths = bam.lengths
    
    # Check if the chromosome exists
    if chromosome in references:
        index = references.index(chromosome)
        chrom_length = lengths[index]
    else:
        chrom_length = None
    
    bam.close()
    return chrom_length

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

def count_stats_for_block(args):
	"""
	Counts statistics for a given block in the chromosome.

	Args:
		args (tuple): Tuple containing (bam_file, chrom, start, end).

	Returns:
		dict: Dictionary with statistics for the block.
	"""
	bam_file, chrom, start, end, outfile = args

	# Open BAM file
	bam = pysam.AlignmentFile(bam_file, "rb")

	unmapped_mate = 0
	mate_diff_chr = 0
	hard_clipped = 0
	None_seq = 0
	total_reads = 0

	oer_reads = []

	dc_r1 = {}
	dc_r2 = {}
	fw = open(outfile,'w')
	# Fetch reads for the block
	for read in bam.fetch(chrom, start, end):
		total_reads += 1

		# Count reads where the mate is unmapped
		if read.mate_is_unmapped:
			unmapped_mate += 1
			oer_reads.append(read.qname)
			continue

		# Count reads where the mate maps to a different chromosome
		if read.next_reference_id != read.reference_id and not read.mate_is_unmapped:
			mate_diff_chr += 1
			oer_reads.append(read.qname)
			continue

		# Count hard-clipped reads
		if "H" in read.cigarstring:
			hard_clipped += 1
			oer_reads.append(read.qname)
			continue

		# read seq is None:
		if read.seq is None:
			None_seq+=1
			oer_reads.append(read.qname)
			continue
		qname = read.qname

		if read.is_reverse:
			seq = reverse_complement(read.seq)
			qual = read.qual[::-1]
		else:
			seq = read.seq
			qual = read.qual
		found_pair = 0
		if read.is_read1:
			if qname in dc_r2:
				seq1, qual1 = seq, qual
				seq2, qual2 = dc_r2[qname]
				dc_r2.pop(qname)
				found_pair = 1
			else:
				dc_r1[qname] = (seq, qual)
		else:
			if qname in dc_r1:
				seq2, qual2 = seq, qual
				seq1, qual1 = dc_r1[qname]
				dc_r1.pop(qname)
				found_pair = 1
			else:
				dc_r2[qname] = (seq, qual)
		if found_pair:
			assert len(seq1) == len(seq2) == len(qual1) == len(qual2)
			fw.write("@"+qname+'\n')
			fw.write(seq1+"\n")
			fw.write("+\n")
			fw.write(qual1+"\n")
			fw.write("@"+qname+'\n')
			fw.write(seq2+"\n")
			fw.write("+\n")
			fw.write(qual2+"\n")
	bam.close()
	fw.close()

	return ({
		"block": f"{chrom}:{start}-{end}",
		"unmapped_mate": unmapped_mate,
		"mate_diff_chr": mate_diff_chr,
		"hard_clipped": hard_clipped,
		"total_reads": total_reads,
		"none_seq": None_seq
	},list(set(oer_reads)))

def merge_files(files, outfile):
	cmd = "cat "+" ".join(files)+" > "+outfile
	print(cmd)
	Popen(cmd, shell= True).wait()


def split_chromosome(chrom_length, num_threads):
    """
    Splits a chromosome into equal blocks based on its length and the number of threads.

    Args:
        chrom_length (int): Length of the chromosome.
        num_threads (int): Number of threads (blocks) to divide the chromosome into.

    Returns:
        list of tuples: List of (start, end) positions for each block.
    """
    block_size = chrom_length // num_threads
    blocks = []

    for i in range(num_threads):
        start = i * block_size
        end = start + block_size - 1 if i < num_threads - 1 else chrom_length - 1
        blocks.append((start, end))

    return blocks


def gen_stats(bam_file, chromosome, out_dir, prefix, num_threads = 50):
	# Input BAM file and chromosome information
	# bam_file = "/lio/lfs/maiziezhou_lab/maiziezhou_lab/Datasets/CG_datasets/V350159389/align_by_longranger/L01/possorted_bam.bam"  # Path to your BAM file
	# chromosome = "chr22"     # Chromosome name
	# num_threads = 60          # Number of parallel threads
	# out_dir = "V350159389_chr22"
	# prefix = "V350159389"
      
	temp_dir = out_dir+"/temp/"
	outfile = f"{out_dir}/{prefix}_{chromosome}_temp.fq"


	chrom_length = get_chromosome_length(bam_file, chromosome)

	os.system("mkdir -p "+ temp_dir)

	# Split chromosome into blocks
	blocks = split_chromosome(chrom_length, num_threads)

	temp_files = [temp_dir + f"/{chromosome}_{start}_{end}.fq" for start, end in blocks]

	# Prepare arguments for parallel processing
	tasks = [(bam_file, chromosome, start, end, temp_dir + f"/{chromosome}_{start}_{end}.fq") for start, end in blocks]
	

	# Process blocks in parallel
	results = Parallel(n_jobs=num_threads)(delayed(count_stats_for_block)(task) for task in tqdm(tasks))
	merge_files(temp_files, outfile)


	# Aggregate results
	total_unmapped_mate = 0
	total_mate_diff_chr = 0
	total_hard_clipped = 0
	total_none_seq = 0
	total_reads = 0

	print("Block-level statistics:")
	all_oer = []
	for result,oer in results:
		print(result)
		total_unmapped_mate += result["unmapped_mate"]
		total_mate_diff_chr += result["mate_diff_chr"]
		total_hard_clipped += result["hard_clipped"]
		total_none_seq += result["none_seq"]
		total_reads += result["total_reads"]
		all_oer.extend(oer)

	# Print summary statistics
	print("\nSummary statistics:")
	print(f"Total unmapped mate: {total_unmapped_mate}")
	print(f"Total mate diff chromosome: {total_mate_diff_chr}")
	print(f"Total hard clipped: {total_hard_clipped}")
	print(f"Total none seq: {total_none_seq}")
	print(f"Total reads: {total_reads}")
	perct = round((total_unmapped_mate + total_mate_diff_chr + total_hard_clipped + total_none_seq)/total_reads * 100,2)
	print(f"Percentage of OER: {perct}%")

	# Open the log file in append mode
	log_file = out_dir+"/summary_statistics.log"
	with open(log_file, 'a') as log:
		# Write the summary statistics
		log.write("\nSummary statistics:\n")
		log.write(f"Total unmapped mate: {total_unmapped_mate}\n")
		log.write(f"Total mate diff chromosome: {total_mate_diff_chr}\n")
		log.write(f"Total hard clipped: {total_hard_clipped}\n")
		log.write(f"Total none seq: {total_none_seq}\n")
		log.write(f"Total reads: {total_reads}\n")
		log.write(f"Percentage of OER: {perct}%\n")


	with open(out_dir+"/oer_"+chromosome,'w' ) as f:
		f.write('\n'.join(list(set(all_oer)))+'\n')

def read_qname_list(qname_file):
    """
    Reads the QNAME list from a file.

    Args:
        qname_file (str): Path to the QNAME file.

    Returns:
        set: Set of read names to extract.
    """
    with open(qname_file, "r") as f:
        return set(line.strip() for line in f)


def process_block(args):
    """
    Extracts reads from a block of the FASTQ file and writes to a temporary file.

    Args:
        args (tuple): Tuple containing (fastq_file, start, end, qname_set, temp_dir).

    Returns:
        str: Path to the temporary file containing matched reads.
    """
    fastq_file, start, end, qname_set, temp_dir = args
    temp_file = tempfile.NamedTemporaryFile(delete=False, dir=temp_dir, suffix=".fastq")
    temp_file_path = temp_file.name
    temp_file.close()  # Close for writing using gzip

    open_func = gzip.open if fastq_file.endswith(".gz") else open
    with open_func(fastq_file, "rt") as f, open(temp_file_path, "w") as out_f:
        # Seek to the start position
        f.seek(start)

        # Adjust position to the next complete FASTQ record
        if start != 0:
            f.readline()  # Skip to the next line


        while f.tell() < end:
            # Read one FASTQ record (4 lines)
            header = f.readline().strip()
            if not header or not header.startswith('@'):
                continue  # Stop if the line doesn't start with '@' (invalid header)
            sequence = f.readline().strip()
            plus = f.readline().strip()
            quality = f.readline().strip()

            if not header or not sequence or not plus or not quality:
                break

            # Check if the read name is in the QNAME set
            read_name = header.split()[0].lstrip("@").split('#')[0]
            if read_name in qname_set:
                out_f.write(f"{header}\n{sequence}\n{plus}\n{quality}\n")

    return temp_file_path


def split_file_by_bytes(fastq_file, num_blocks):
    """
    Splits a FASTQ file into byte blocks for parallel processing.

    Args:
        fastq_file (str): Path to the FASTQ file.
        num_blocks (int): Number of blocks to divide the file into.

    Returns:
        list of tuples: List of (start_byte, end_byte) for each block.
    """
    file_size = os.path.getsize(fastq_file)
    block_size = file_size // num_blocks

    blocks = []
    for i in range(num_blocks):
        start = i * block_size
        end = (start + block_size - 1) if i < num_blocks - 1 else file_size - 1
        blocks.append((start, end))

    return blocks


def merge_temp_files(temp_files, output_file):
    """
    Merges temporary files into a single output file.

    Args:
        temp_files (list of str): List of paths to temporary files.
        output_file (str): Path to the output FASTQ file.
    """
    cmd = "cat " + " ".join(temp_files) + " > " + output_file
    print(cmd)
    Popen(cmd, shell= True).wait()
    


def extract_oer(fastq_file, out_dir, chromosome, prefix, num_threads = 50):
    # Input files and parameters
    # fastq_file = "/lio/lfs/maiziezhou_lab/maiziezhou_lab/Datasets/CG_datasets/V350159389/align_by_longranger/L01/test_sample_S1_L001_R1_001.fastq"  # Path to the WGS FASTQ file
    # out_dir = "V350159389_chr22"
    # prefix = "V350159389"
    # chromosome = "chr22"
    # num_threads = 50  # Number of threads
    
	outfile = f"{out_dir}/{prefix}_{chromosome}_temp.fq"

	qname_file = out_dir+"/oer_"+chromosome
	oer_output_file = qname_file + ".fq"  # Output FASTQ file for extracted reads


	# Create a temporary directory
	temp_dir = out_dir+"/temp/"
	# temp_dir = "./temp_dir/"
	print(temp_dir)

	try:
		# Read QNAME list
		qname_set = read_qname_list(qname_file)

		# Split the FASTQ file into blocks
		blocks = split_file_by_bytes(fastq_file, num_threads)

		# Prepare arguments for parallel processing
		tasks = [(fastq_file, start, end, qname_set, temp_dir) for start, end in blocks]

		# Process blocks in parallel
		with Pool(num_threads) as pool:
			temp_files = pool.map(process_block, tasks)

		# Merge temporary files into the final output
		merge_temp_files(temp_files, oer_output_file)

		print(f"Extraction complete. Output written to {oer_output_file}")

		cmd = f"cat {oer_output_file} >> {outfile}"
		Popen(cmd, shell= True).wait()
		

	finally:
		# Clean up temporary directory
		os.system("rm -r "+temp_dir)

def remove_duplicates(input_fastq, output_fastq, out_dir):
	"""
	Remove duplicate read pairs from a FASTQ file.

	Parameters:
		input_fastq (str): Path to the input FASTQ file.
		output_fastq (str): Path to the output FASTQ file without duplicates.
	"""
	seen = set()  # To store unique read names
	total_pair = 0
	dup_pair = 0
	cnt = 0 
	with open(input_fastq, 'r') as infile, open(output_fastq, 'w') as outfile:
		lines = []  # Buffer for the current read pair
		for i, line in enumerate(infile):
			lines.append(line)
			if i  % 8 == 7:  # Every 8 lines represent a read pair
				total_pair+=1
				read_name = lines[0].split()[0]  # Extract read name from the first line
				if read_name not in seen:
					seen.add(read_name)
					outfile.writelines(lines)  # Write the unique read pair to output
					cnt+=1
				else:
					dup_pair+=1
				lines = []  # Clear buffer for the next read pair
	pct = round(dup_pair/total_pair*100,2)
	with open(out_dir + "/dedup_summary.txt",'w') as f:
		f.write(f"Total pair of reads: {total_pair}\n" )
		f.write(f"Duplicated pair of reads: {dup_pair}({pct}%)\n" )
		f.write(f"final dedupped pair of reads: {cnt}\n" )

		

def main(bam_file, fastq_file, chromosome,out_dir, prefix,num_threads):
	gen_stats(bam_file, chromosome, out_dir, prefix, num_threads )
	extract_oer(fastq_file, out_dir, chromosome,  prefix, num_threads )
	input_fastq =  f"{out_dir}/{prefix}_{chromosome}_temp.fq"
	output_fastq =  f"{out_dir}/{prefix}_{chromosome}.fq"
	remove_duplicates(input_fastq, output_fastq, out_dir)



def parse_args():
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(description="Process BAM and FASTQ files with specific parameters.")
    
    parser.add_argument(
        "--bam_file", '-bam',
        type=str, 
        required=True, 
        help="Path to the input BAM file."
    )
    parser.add_argument(
        "--fastq_file", '-fq',
        type=str, 
        required=True, 
        help="Path to the WGS FASTQ file."
    )
    parser.add_argument(
        "--chromosome", '-chr',
        type=str, 
        required=True, 
        help="Chromosome name (e.g., chr6)."
    )
    parser.add_argument(
        "--num_threads", '-t',
        type=int, 
        default=30, 
        help="Number of parallel threads to use. Default is 30."
    )
    parser.add_argument(
        "--out_dir", '-o',
        type=str, 
        required=True, 
        help="Output directory for results."
    )
    parser.add_argument(
        "--prefix", '-px',
        type=str, 
        required=True, 
        help="Prefix for output files."
    )
    
    return parser.parse_args()

if __name__ == "__main__":
	args = parse_args()
	# Print parsed arguments
	print(f"BAM file: {args.bam_file}")
	print(f"FASTQ file: {args.fastq_file}")
	print(f"Chromosome: {args.chromosome}")
	print(f"Number of threads: {args.num_threads}")
	print(f"Output directory: {args.out_dir}")
	print(f"Prefix: {args.prefix}")



	# bam_file = "/lio/lfs/maiziezhou_lab/maiziezhou_lab/Datasets/CG_datasets/V350159389/align_by_longranger/L01/possorted_bam.bam"  # Path to your BAM file
	# fastq_file = "/lio/lfs/maiziezhou_lab/maiziezhou_lab/Datasets/CG_datasets/V350159389/align_by_longranger/L01/test_sample_S1_L001_R1_001.fastq"  # Path to the WGS FASTQ file
	# chromosome = "chr6"     # Chromosome name
	# num_threads = 60          # Number of parallel threads
	# out_dir = "V350159389_L01_chr6"
	# prefix = "V350159389"
	main(args.bam_file, args.fastq_file, args.chromosome,args.out_dir, args.prefix,args.num_threads)
