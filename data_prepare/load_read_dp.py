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
    # print(st)
    # print(file)
    while line:=f.readline():
        # print(line)
        cur_rnv,tst,ted = line.split()
        if cur_rnv == rnv:
            return (int(tst), int(ted))
        if f.tell()>= ed:
            break 

    print(f"could not file {rnv} in {file}")
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







def trace_rnv_in_tree(indir,n_idx,rnv, idxl):
    st  = 0
    ed = os.path.getsize(indir+'/index_%d.txt'%(n_idx-1))
    reuse_list = []
    for i in range(n_idx):
        k = n_idx - i - 1
        infile = indir+'/index_%d.txt'% k
        st,ed = locate_rnv(infile, rnv[:i+1],st,ed )
        # print(st,ed,indir+'/index_%d.txt'% (k-1))
        reuse_list.append((st,ed))

    rg_list = locate_merged_index(indir,rnv[:-1],st,ed)
    # print(rg_list)

    reuse_list.append(rg_list)

    result = locate_final_index(indir, rnv, rg_list, idxl)
    reuse_list.append(result)
    # print(result)
    return result, reuse_list


def trace_rnv_in_tree_with_prior(indir,n_idx,rnv, idxl, reuse_list,n_equal):
    if n_equal == 0:
        st  = 0
        ed = os.path.getsize(indir+'/index_%d.txt'%(n_idx-1))
        reuse_list = []
        for i in range(n_idx):
            k = n_idx - i - 1
            infile = indir+'/index_%d.txt'% k
            st,ed = locate_rnv(infile, rnv[:i+1],st,ed )
            reuse_list.append((st,ed))

        rg_list = locate_merged_index(indir,rnv[:-1],st,ed)
        # print(rg_list)

        reuse_list.append(rg_list)

        result = locate_final_index(indir, rnv, rg_list, idxl)
        reuse_list.append(result)
        # print(result)
        return result, reuse_list
    elif n_equal == (len(idxl)-1):
        rg_list = reuse_list[-2]
        result = locate_final_index(indir, rnv, rg_list, idxl)
        reuse_list[-1] = result
        # print(result)
        return result, reuse_list
    elif n_equal == (len(idxl)-2):
        st,ed = reuse_list[-3]

        rg_list = locate_merged_index(indir,rnv[:-1],st,ed)
        # print(rg_list)

        reuse_list[-2] = rg_list

        result = locate_final_index(indir, rnv, rg_list, idxl)
        reuse_list[-1] = result
        return result, reuse_list
    else:
        # print('her')
        # print(reuse_list)
        # print(n_equal)
        st,ed  = reuse_list[n_equal-1]
        # ed = os.path.getsize(indir+'/index_%d.txt'%(n_idx-1))
        reuse_list = reuse_list[:n_equal].copy()
        # print(len(reuse_list))
        for i in range(n_equal,n_idx):
            k = n_idx - i - 1
            # print(k)
            infile = indir+'/index_%d.txt'% k
            # print(infile,st,ed)
            st,ed = locate_rnv(infile, rnv[:i+1],st,ed )
            reuse_list.append((st,ed))

        rg_list = locate_merged_index(indir,rnv[:-1],st,ed)
        # print(rg_list)

        reuse_list.append(rg_list)

        result = locate_final_index(indir, rnv, rg_list, idxl)
        reuse_list.append(result)
        # print(result)
        return result, reuse_list






def extract_read(tfile,st,ed):
    f = open(tfile,'r')
    f.seek(st)
    lines = []
    while line:=f.readline():
        lines.append(line)
        if f.tell()>= ed:
            break 
    return ''.join(lines)

def convert_rn(read_name,idxl):
    return ''.join([read_name[i] for i in idxl ])

def retrive_read(read_name,idxl, indir, tfile):
    target_rnv  = ''.join([read_name[i] for i in idxl ])
    # print(target_rnv)
    n_idx = len(idxl) -2 
    result, reuse_list = trace_rnv_in_tree(indir,n_idx,target_rnv, idxl)
    st,ed = result
    lines = extract_read(tfile,st,ed)
    return lines, reuse_list

def retrive_read_with_prior(read_name,idxl, indir, tfile,reuse_list, n_equal):
    target_rnv  = ''.join([read_name[i] for i in idxl ])
    # print(target_rnv)
    n_idx = len(idxl) -2 
    result, reuse_list = trace_rnv_in_tree_with_prior(indir,n_idx,target_rnv, idxl,reuse_list, n_equal)
    st,ed = result
    lines = extract_read(tfile,st,ed)
    return lines, reuse_list

def check_qual(rn1,rn2):
    cnt =0
    for i in range(len(rn1)):
        if rn1[i]==rn2[i]:
            cnt+=1
        else:
            break 

    return cnt



def retrieve_reads(rn_list,indir,tfile, outfile):
    fw = open(outfile,'w')
    var_file = indir +'/var_index.txt'
    idxl = load_var(var_file)

    read, reuse_list = retrive_read(rn_list[0],idxl, indir, tfile) 
    fw.write(read)

    for i in tqdm(range(1,len(rn_list))):
        rn1 = convert_rn(rn_list[i-1],idxl)
        rn2 = convert_rn(rn_list[i],idxl)
        n_equal = check_qual(rn1,rn2)
        # print(rn_list[i],'n e',n_equal)
        read, reuse_list = retrive_read_with_prior(rn_list[i],idxl, indir, tfile,reuse_list, n_equal)
        fw.write(read)

    return 



# # read_name = "V350158968L1C003R03800068259"
# tfile = "/lio/lfs/maiziezhou_lab/maiziezhou_lab/Datasets/CG_datasets/V350158968/merged_split_read.1.fq"
# indir = "/data/maiziezhou_lab/CanLuo/Software/Fastq_Tool/out_r1/"
# outfile = 'test.fq'

# names = []
# i = 0
# with open("/data/maiziezhou_lab/ZhouLab_Projects/VandyCG_stLFR/Experiments/2023_October/Aquila_stLFR_chr6/read/unused_name_uniq_sorted.txt",'r') as f:
#     for line in f:
#         i+=1
#         if i>1000:
#             break
#         names.append(line[:-1])


# retrieve_reads(names,indir,tfile, outfile)





# lines = retrive_read(read_name,idxl, indir, tfile)











































