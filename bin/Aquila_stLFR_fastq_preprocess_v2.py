#!/usr/bin/env python
import gzip
from argparse import ArgumentParser
import sys
parser = ArgumentParser(description="Preprocessing paired fastq files for Aquila_stLFR\n",usage='use "python3 %(prog)s --help" for more information')
parser.add_argument('--fastq_1','-1',help="origin stLFR fastq 1 (gz file)",required=True)
parser.add_argument('--fastq_2','-2',help="origin stLFR fastq 2 (gz file)",required=True)
parser.add_argument('--out_file','-o',help="output stLFR fastq file for Aquila_stLFR",default="stLFR_fastq_for_Aquila.fastq")
args = parser.parse_args()


def merge_paired_reads_gz(read1_file,read2_file,out_file):
    fw = open(out_file,"w")
    count = 0
    print("v2")
    def smart_open(filename, mode="rt"):
        """Open a file normally if uncompressed, otherwise use gzip.open"""
        return gzip.open(filename, mode) if filename.endswith(".gz") else open(filename, mode)

    with smart_open(read1_file,"r") as f1, smart_open(read2_file,"r") as f2:
        for line1,line2 in zip(f1,f2):
            if count%4 == 0:
                write_flag = 1
                data1 = line1.decode().rsplit()
                data2 = line2.decode().rsplit()
                qname1,bx = data1[0].split("#")
                print(data1,qname1,bx)
                exit()
                # bx = data1[1]
                # bc = data1[2]
                if bx == '0_0_0':
                    write_flag = 0
                qname2 = data2[0].split("#")[0]
                if write_flag:
                    fw.writelines(qname1 + "\tBX:Z:" + bx + "\n")
                line5 = qname2 + "\tBX:Z" + bx +  "\n"
            elif count%4 == 1:
                if write_flag:
                    fw.writelines(line1.decode())
                line6 = line2.decode() 
            elif count%4 == 2:
                if write_flag:
                    fw.writelines(line1.decode())
                line7 = line2.decode() 
            elif count%4 == 3:
                line8 = line2.decode() 
                if write_flag:
                    fw.writelines(line1.decode())  
                    fw.writelines(line5)
                    fw.writelines(line6)
                    fw.writelines(line7)
                    fw.writelines(line8)
            count += 1                
    fw.close()


def merge_paired_reads(read1_file,read2_file,out_file):
    fw = open(out_file,"w")
    count = 0

    def smart_open(filename, mode="rt"):
        """Open a file normally if uncompressed, otherwise use gzip.open"""
        return gzip.open(filename, mode) if filename.endswith(".gz") else open(filename, mode)

    with smart_open(read1_file,"r") as f1, smart_open(read2_file,"r") as f2:
        for line1,line2 in zip(f1,f2):
            if count%4 == 0:
                write_flag = 1
                data1 = line1.rsplit()
                data2 = line2.rsplit()
                qname1,bx = data1[0].split("#")
                bx = bx.split("/")[0]
                # bx = data1[1]
                # bc = data1[2]
                if bx == '0_0_0':
                    write_flag = 0
                qname2 = data2[0].split("#")[0]
                if write_flag:
                    fw.writelines(qname1 + "\tBX:Z:" + bx + "\n")
                line5 = qname2 + "\tBX:Z:" + bx  + "\n"
            elif count%4 == 1:
                if write_flag:
                    fw.writelines(line1)
                line6 = line2
            elif count%4 == 2:
                if write_flag:
                    fw.writelines(line1)
                line7 = line2 
            elif count%4 == 3:
                line8 = line2
                if write_flag:
                    fw.writelines(line1)  
                    fw.writelines(line5)
                    fw.writelines(line6)
                    fw.writelines(line7)
                    fw.writelines(line8)
            count += 1                
    fw.close()



def main():
    if len(sys.argv) == 1:
        Popen("python3 " + "Aquila_stLFR_fastq_preprocess.py -h",shell=True).wait()
    else:
        fastq_1 = args.fastq_1    
        fastq_2 = args.fastq_2    
        out_file = args.out_file  
        if fastq_1.endswith(".gz"):
            merge_paired_reads_gz(fastq_1,fastq_2,out_file)
        else:
            merge_paired_reads(fastq_1,fastq_2,out_file)


if __name__ == "__main__":
    main()

