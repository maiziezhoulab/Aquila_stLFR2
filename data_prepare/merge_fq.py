import argparse

# def merge_fastq(input1, input2, output):
#     with open(input1, 'r') as file1, open(input2, 'r') as file2, open(output, 'w') as output_file:
#         for line1, line2 in zip(file1, file2):
#             output_file.write(line1)
#             output_file.write(line2)

def merge_fastq(input1, input2, output):
    with open(input1, 'r') as file1, open(input2, 'r') as file2, open(output, 'w') as output_file:
        while True:
            try:
                for _ in range(4):
                    if _ == 0:
                        # line = 
                        output_file.write(next(file1).split('#')[0]+'\n')
                    else:

                        output_file.write(next(file1))
                for _ in range(4):
                    if _ == 0:
                         output_file.write(next(file2).split('#')[0]+'\n')
                    else:
                        output_file.write(next(file2))
            except StopIteration:
                break

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Merge two paired FASTQ files into a single interleaved file.')
    parser.add_argument('--input1','-i1', help='Input FASTQ file 1')
    parser.add_argument('--input2', '-i2',help='Input FASTQ file 2')
    parser.add_argument('--output', '-o',help='Output interleaved FASTQ file')

    args = parser.parse_args()

    merge_fastq(args.input1, args.input2, args.output)
