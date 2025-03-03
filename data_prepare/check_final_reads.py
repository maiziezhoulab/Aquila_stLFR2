def load_read_names_from_fastq(fastq_file):
    read_names = []
    i = 0
    with open(fastq_file, 'r') as f:
        for line in f:
            i+=1
            if i%4 == 1:
                # Extract the read name from the FASTQ header line
                read_name = line.strip().split()[0][1:]
                read_names.append(read_name)
                if i%8 ==  5 :
                    assert read_name ==  read_names[-1]
                    
    return set(read_names)

def load_read_names_from_file(read_name_file):
    read_names = []
    with open(read_name_file, 'r') as f:
        for line in f:
            read_names.append(line.strip())
    return set(read_names)

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Load read names from a FASTQ file and a read name file.")
    parser.add_argument("--fastq_file",'-fq', help="FASTQ file to extract read names from")
    parser.add_argument("--read_name_file",'-rn', help="File containing read names (one per line)")

    args = parser.parse_args()

    fastq_read_names = load_read_names_from_fastq(args.fastq_file)
    read_name_file_read_names = load_read_names_from_file(args.read_name_file)

    print("uniq names in fq:",len(fastq_read_names))
    print("uniq names in rn:",len(read_name_file_read_names))

    diff_read_names = read_name_file_read_names - fastq_read_names
    print("in rn but not in fq:", len(diff_read_names))
    assert len(diff_read_names) == 0

