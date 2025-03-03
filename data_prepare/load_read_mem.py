import os
import pandas as pd
# from joblib import Parallel, delayed
from tqdm import tqdm
import numpy as np
from collections import defaultdict
from subprocess import Popen
def load_var(var_file):

    with open(var_file,'r') as f:
        idxl = eval(f.read())

    return idxl 


def locate_rnv(file, rnv, st,ed):
    # f = open(file,'r')
    f = file
    f.seek(st)

    while line:=f.readline():
        cur_rnv,tst,ted = line.split()
        if cur_rnv == rnv:
            return (int(tst), int(ted))
        if f.tell()>= ed:
            break 

    print(f"could not find {rnv} in {file}")
    exit()


def locate_merged_index(f, rnv,st,ed):
    # infile = f"{indir}/merged_first_index_merge_entry.txt"

    # f = open(infile,'r')
    f.seek(st)

    while line:=f.readline():
        cur_rnv = line.split()[0]
        if cur_rnv == rnv:
            return [int(x) for x in line.split()[1:]]
        if f.tell()>= ed:
            break 
    print(f"could not find {rnv} in merged_first_index_merge_entry.txt ")
    exit()

def extract_rnv(rn,idxl):
    return ''.join([rn[i] for i in idxl])

def locate_final_index_one_file(f,rnv,st,ed,idxl):

    # f = open(infile,'r')
    f.seek(st)
    while line:=f.readline():
        cur_rnv = extract_rnv(line.split()[0], idxl)
        if cur_rnv == rnv:
            return [int(x) for x in line.split()[1:]]
        if f.tell()>= ed:
            break 
    return -1


def locate_final_index(f_list, rnv, rg_list, idxl):
    # df = pd.read_csv(indir+'/chunks.csv')

    for i in range(0,len(rg_list),3):
        fid, st,ed = rg_list[i:i+3]
        # fst = df['start'][fid]
        # fed = df['end'][fid]
        # infile = f'{indir}/NA24385_{fst}_{fed}_sorted.txt'
        infile = f_list[fid]
        result = locate_final_index_one_file(infile,rnv,st,ed,idxl)
        if result!= -1:
            return result
    print(f"could not find {rnv} in final index from ",rg_list)
    exit()


def open_all_final_index(indir, out_prefix, ):
    df = pd.read_csv(indir+'/../chunks.csv')
    f_list = []

    for fid in range(df.shape[0]):
        fst = df['start'][fid]
        fed = df['end'][fid]
        infile = f'{indir}/{out_prefix}_{fst}_{fed}_sorted.txt'
        f = open(infile,'r')
        f_list.append(f)

    return f_list
        






def trace_rnv_in_tree(indir,n_idx,rnv, idxl,st,ed,f_id1, f_id0,f_mg, f_list):
    # st  = 0
    # ed = os.path.getsize(indir+'/index_%d.txt'%(n_idx-1))
    # print(n_idx)
    for i in range(n_idx - 2,n_idx):
        k = n_idx - i - 1
        if k == 1:
            infile =f_id1
        elif k == 0:
            infile = f_id0
        # print(infile, st, ed)
        st,ed = locate_rnv(infile, rnv[:i+1],st,ed )


    rg_list = locate_merged_index(f_mg,rnv[:-1],st,ed)
    # print(rg_list)

    result = locate_final_index(f_list, rnv, rg_list, idxl)
    # print(result)
    return result


def extract_read(tfile,rg_list,outfile):
    with open(outfile,'w') as fw:
        f = open(tfile,'r')
        for rg in tqdm(rg_list, desc = f"write to chunk {outfile}"):
            st = rg[0]
            f.seek(st)
            for i in range(4):
                l = f.readline()
                fw.write(l)
        f.close()

    # lines = []
    # while line:=f.readline():
    #     lines.append(line)
    #     if f.tell()>= ed:
    #         break 
    # return ''.join(lines)
    # return -1


def retrive_read(read_name,idxl, indir, tfile, dc_idx2,f_id1, f_id0,f_mg, f_list):
    target_rnv  = ''.join([read_name[i] for i in idxl ])
    # print(target_rnv)
    n_idx = len(idxl) -2 
    st,ed = dc_idx2[target_rnv[:-4]]
    st,ed = trace_rnv_in_tree(indir,n_idx,target_rnv, idxl, st,ed,f_id1, f_id0,f_mg, f_list)
    # lines = extract_read(tfile,st,ed)
    return st,ed


def load_idx2(indx2):
    dc = {}

    with open(indx2,'r')as f:
        for line in f:
            data = line.split()
            dc[data[0]] = (int(data[1]), int(data[2]))

    return dc 

def retrieve_reads(rn_list,indir,tfile, outfile, out_prefix = "Sample"):
    var_file = indir +'/var_index.txt'
    idxl = load_var(var_file)
    dc_idx2 = load_idx2(indir+'/index_2.txt')
    # print(len(dc_idx2))
    rg_list = []

    #--------- load file first
    # print("open all index...")
    f_id1 = open(indir + '/index_1.txt')
    f_id0 = open(indir + '/index_0.txt')
    f_mg = open(indir + '/merged_first_index_merge_entry.txt')
    f_list = open_all_final_index(indir, out_prefix)
    # print("finish opening all index")


    # with open(outfile,'w') as fw:
    for read_name in tqdm(rn_list):
        rg = retrive_read(read_name,idxl, indir, tfile, dc_idx2, f_id1, f_id0,f_mg, f_list) 
        rg_list.append(rg)
            # for rg in rg_list:
            # fw.write(str(rg[0])+','+str(rg[1])+'\n')
            # fw.write(read)
    extract_read(tfile, rg_list, outfile)

    #--------- close all index
    # print("close all index...")
    f_id1.close()
    f_id0.close()
    f_mg.close()
    for f in f_list:
        f.close()
    # print("finish closing all index")
    return 


def split_one_list(inlist):
    dc = defaultdict(list)
 
    for rn in inlist:
        k = len(rn)
        dc[k].append(rn)
    return dc
        


def retrieve_reads_mt_len(inlist, indir, tfile, outfile, out_prefix, temp_dir,temp_prefix):
    # temp_dir = outdir+'/temp/'
    if not os.path.exists(temp_dir):
        os.system("mkdir -p "+ temp_dir)

    dc_k = split_one_list(inlist)
    temp_fq_files = [temp_dir+f'/{temp_prefix}_rn_len{k}.fq'   for k in dc_k]
    k_list = list(dc_k.keys())
    sub_dir_list = [ indir +f'/rn_len{k}/'  for k in k_list]

    #-------------------retrieve reads by qname length

    for i in range(len(k_list)):
        rn_list = dc_k[k_list[i]]
        sub_dir = sub_dir_list[i]
        temp_fq = temp_fq_files[i]
        retrieve_reads(rn_list,sub_dir,tfile, temp_fq , out_prefix)

    #-------------------reduce files
    cmd = f"cp {temp_fq_files[0]} {outfile}"
    Popen(cmd, shell = True).wait()

    for i in range(1,len(k_list)):
        cmd = f"cat {temp_fq_files[i]} >> {outfile}"
        Popen(cmd, shell = True).wait()

    #-------------clean temp files
    cmd = "rm " + ' '.join(temp_fq_files)
    Popen(cmd, shell = True).wait()








    

if __name__ == "__main__":
    import argparse
    from argparse import ArgumentParser
    parser = ArgumentParser(description="",
        usage='use "python3 %(prog)s --help" for more information',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--input_path','-i')
    parser.add_argument('--fastq_path','-fq')
    parser.add_argument('--index_dir','-idx')
    parser.add_argument('--temp_dir','-tmp')
    parser.add_argument('--temp_prefix','-tpx')
    parser.add_argument('--out_prefix','-px', default = "Sample")
    parser.add_argument('--output_path','-o')
    # parser.add_argument('--n_thread','-t', type = int, default = 22 )
    # parser.add_argument('--delete_temp_file','-d', action='store_true')
    args = parser.parse_args()
    input_path = args.input_path
    fastq_path = args.fastq_path
    index_dir = args.index_dir
    temp_dir = args.temp_dir
    temp_prefix = args.temp_prefix
    out_prefix = args.out_prefix

    output_path = args.output_path
    # n_thread = args.n_thread

    with open(input_path,'r') as f:
        rn_list = f.read().split('\n')[:-1]
    print("num of read names:",len(rn_list))

    retrieve_reads_mt_len(rn_list, index_dir,fastq_path, output_path, 
    out_prefix, temp_dir,temp_prefix)


    # retrieve_reads(rn_list,index_dir,fastq_path, output_path, out_prefix = "Sample")










































