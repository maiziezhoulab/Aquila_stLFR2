import os
import pandas as pd
from joblib import Parallel, delayed
from tqdm import tqdm
def load_var(var_file):

    with open(var_file,'r') as f:
        idxl = eval(f.read())

    return idxl 


def locate_rnv(file, rnv, st,ed):
    f = open(file,'r')
    f.seek(st)

    while line:=f.readline():
        cur_rnv,tst,ted = line.split()
        if cur_rnv == rnv:
            return (int(tst), int(ted))
        if f.tell()>= ed:
            break 

    print(f"could not find {rnv} in {file}")
    exit()


def locate_merged_index(indir, rnv,st,ed):
    infile = f"{indir}/merged_first_index_merge_entry.txt"

    f = open(infile,'r')
    f.seek(st)

    while line:=f.readline():
        cur_rnv = line.split()[0]
        if cur_rnv == rnv:
            return [int(x) for x in line.split()[1:]]
        if f.tell()>= ed:
            break 
    print(f"could not file {rnv} in {infile}")
    exit()

def extract_rnv(rn,idxl):
    return ''.join([rn[i] for i in idxl])

def locate_final_index_one_file(infile,rnv,st,ed,idxl):

    f = open(infile,'r')
    f.seek(st)
    while line:=f.readline():
        cur_rnv = extract_rnv(line.split()[0], idxl)
        if cur_rnv == rnv:
            return [int(x) for x in line.split()[1:]]
        if f.tell()>= ed:
            break 
    return -1


def locate_final_index(indir, rnv, rg_list, idxl):
    df = pd.read_csv(indir+'/chunks.csv')

    for i in range(0,len(rg_list),3):
        fid, st,ed = rg_list[i:i+3]
        fst = df['start'][fid]
        fed = df['end'][fid]
        infile = f'{indir}/NA24385_{fst}_{fed}_sorted.txt'
        result = locate_final_index_one_file(infile,rnv,st,ed,idxl)
        if result!= -1:
            return result
    print(f"could not find {rnv} in final index from ",rg_list)
    exit()







def trace_rnv_in_tree(indir,n_idx,rnv, idxl,st,ed):
    # st  = 0
    # ed = os.path.getsize(indir+'/index_%d.txt'%(n_idx-1))
    for i in range(7,n_idx):
        k = n_idx - i - 1
        infile = indir+'/index_%d.txt'% k
        st,ed = locate_rnv(infile, rnv[:i+1],st,ed )
        if i==(n_idx -1):
            # print(f"Found {rnv[:i+1]} in {infile}, point to {indir}/merged_first_index_merge_entry.txt:{st}-{ed}")
            pass
        else:

            # print(f"Found {rnv[:i+1]} in {infile}, point to {indir+'/index_%d.txt'% (k-1)}:{st}-{ed}")
            pass

    rg_list = locate_merged_index(indir,rnv[:-1],st,ed)
    # print(rg_list)

    result = locate_final_index(indir, rnv, rg_list, idxl)
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


def retrive_read(read_name,idxl, indir, tfile, dc_idx2):
    target_rnv  = ''.join([read_name[i] for i in idxl ])
    # print(target_rnv)
    n_idx = len(idxl) -2 
    st,ed = dc_idx2[target_rnv[:7]]
    st,ed = trace_rnv_in_tree(indir,n_idx,target_rnv, idxl, st,ed)
    # lines = extract_read(tfile,st,ed)
    return st,ed


def load_idx2(indx2):
    dc = {}

    with open(indx2,'r')as f:
        for line in f:
            data = line.split()
            dc[data[0]] = (int(data[1]), int(data[2]))

    return dc 
def retrieve_reads(rn_list,indir,tfile, outfile):
    var_file = indir +'/var_index.txt'
    idxl = load_var(var_file)
    dc_idx2 = load_idx2(indir+'/index_2.txt')
    print(len(dc_idx2))
    rg_list = []
    # with open(outfile,'w') as fw:
    for read_name in tqdm(rn_list):
        rg = retrive_read(read_name,idxl, indir, tfile, dc_idx2) 
        rg_list.append(rg)
            # for rg in rg_list:
            # fw.write(str(rg[0])+','+str(rg[1])+'\n')
            # fw.write(read)
    extract_read(tfile, rg_list, outfile)
    return 



# read_name = "V350158968L1C003R03800068259"
# tfile = "/lio/lfs/maiziezhou_lab/maiziezhou_lab/Datasets/CG_datasets/V350158968/merged_split_read.1.fq"
# indir = "./out_r1/"







# lines = retrive_read(read_name,idxl, indir, tfile)











































