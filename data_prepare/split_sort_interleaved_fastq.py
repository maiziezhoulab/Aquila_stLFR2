

import argparse
from argparse import ArgumentParser
parser = ArgumentParser(description="",
	usage='use "python3 %(prog)s --help" for more information',
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--input_path','-i')
parser.add_argument('--output_dir','-o')
parser.add_argument('--out_prefix','-px')
# parser.add_argument('--n_thread','-t', type = int, default = 22 )
# parser.add_argument('--delete_temp_file','-d', action='store_true')
args = parser.parse_args()
input_path = args.input_path
output_dir = args.output_dir
out_prefix = args.out_prefix






def split_fastq(infile, out1,out2):
    i = 0
    f1 = open(out1,'w')
    f2 = open(out2,'w')

    f = open(infile,'r')
    filesize = os.path.getsize(infile)
    print("input file size: %d GB" %(filesize/(1024**3)))
    while line:= f.readline():
        if i%10000000==0:
            print(f"processed {i} lines ({round(f.tell()/(1024**3),2)}GB)")
        if i%8 in {0,4}:
            line = line.split()[0].split('#')[0]+'\n'
        if i%8 <4:
            f1.write(line)
        else:
            f2.write(line)

        i+=1

    f.close()
    f1.close()
    f2.close()

    

def sort_one_fq(infq,outfq):
    # cmd = '''cat %s | awk 'ORS=NR%%4?"\t":"\n"' |  sort -t$'\t' -k1,1 | awk -F'\t' '{ print $1; print $2; print $3; print $4}' > %s'''%(infq, outfq)
    cmd = '''cat %s | awk 'ORS=NR%%4?"\\t":"\\n"' | sort -t$'\\t' -k1,1 | awk -F'\\t' '{ print $1; print $2; print $3; print $4}' > %s''' % (infq, outfq)
    print(cmd)
    Popen(cmd, shell= True).wait()


def dedup(infq, outfq):

    i = 0
    fw = open(outfq,'w')


    f = open(infq,'r')
    filesize = os.path.getsize(infq)
    print("input file size: %d GB" %(filesize/(1024**3)))
    last_rn = 'heihei'
    dup_cnt = 0
    uniq_cnt = 0
    while line:= f.readline():
        if i%10000000==0:
            print(f"processed {i} lines ({round(f.tell()/(1024**3),2)}GB)")
        if i%4 ==0:
            cur_rn = line.split()[0].split('#')[0]

            if cur_rn!=last_rn:
                write = 1
                uniq_cnt +=1
            else:
                write = 0
                dup_cnt +=1

            last_rn = cur_rn

        if write:
            fw.write(line)



        i+=1

    f.close()
    fw.close()
    
    print(f"{infq} unqi count= {uniq_cnt} dup count= {dup_cnt} dup rate = {round(dup_cnt/uniq_cnt*100,2)}%%")



import os
from subprocess import Popen
from joblib import Parallel, delayed
from tqdm import tqdm

outfile1 = output_dir + '/'+out_prefix +'_R1.fastq'
outfile2 = output_dir + '/'+out_prefix +'_R2.fastq'
outfile1_st = output_dir + '/'+out_prefix +'_R1_sorted.fastq'
outfile2_st = output_dir + '/'+out_prefix +'_R2_sorted.fastq'

outfile1_dedup = output_dir + '/'+out_prefix +'_R1_sorted_dedup.fastq'
outfile2_dedup = output_dir + '/'+out_prefix +'_R2_sorted_dedup.fastq'
in_list = [outfile1,outfile2 ]
out_list = [outfile1_st,outfile2_st ]
out_list_dedup = [outfile1_dedup,outfile2_dedup]

n_thread = 2


if not os.path.exists(output_dir):
    os.system("mkdir -p "+ output_dir)
    
# split_fastq(input_path, outfile1, outfile2)

# sequences = Parallel(n_jobs=n_thread)(delayed(sort_one_fq)(in_list[i], out_list[i]) for i in range(2))

sequences = Parallel(n_jobs=n_thread)(delayed(dedup)(out_list[i], out_list_dedup[i]) for i in range(2))

