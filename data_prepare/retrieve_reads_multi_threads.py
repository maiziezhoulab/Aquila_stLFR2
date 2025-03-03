
if __name__ == "__main__":
    import argparse
    from argparse import ArgumentParser
    parser = ArgumentParser(description="",
        usage='use "python3 %(prog)s --help" for more information',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--fastq_file','-fq')
    parser.add_argument('--qname_file','-q')
    parser.add_argument('--qname_separator','-qsep', default= None, help = "qname field separator;optional, if not provided, will use anything before the first tab as qname in the qname line")
    parser.add_argument('--output_path','-o')
    parser.add_argument('--index_dir','-i')
    # parser.add_argument('--var_path','-v')
    parser.add_argument('--n_thread','-t', type = int, default = 22 )
    parser.add_argument('--delete_temp_file','-d', action='store_true')
    args = parser.parse_args()
    fastq_file = args.fastq_file
    qname_file = args.qname_file
    output_path = args.output_path
    index_dir = args.index_dir
    qname_separator = args.qname_separator
    # var_path = args.var_path
    n_thread = args.n_thread




import pandas as pd
from tqdm import tqdm
from collections import defaultdict
import pickle
import numpy as np
from joblib import Parallel, delayed
from subprocess import Popen
def multi_all(inlist):
    x = 1
    for y in inlist:
        x = x*y
    return x


def renumber_tree(tree):
    new_tree = []
    for x in tree:
        new_tree.append(list(range(len(x))))
    return new_tree

def get_transformer_level1(tree, new_tree):
    trans_list = []
    for i in range(len(tree)):
        trans_list.append(dict(zip(tree[i],new_tree[i])))
    return trans_list

def print_tree(tree):
    print("\n--------------------")
    for i in range(len(tree)):
        print(i,tree[i])
    print("--------------------\n")
    return 

def load_var(infile):
    var_index_list = []
    var_len_list = []
    var_list = []
    with open(infile, 'r') as f:
        for line in f:
            data = line.split()
            if len(data) > 2:
                var_index_list.append(int(data[0]))
                var_len_list.append(len(data)-1)
                var_list.append(data[1:])

    

    size_list = [ multi_all(var_len_list[i+1:])   for i in range(len(var_len_list)-1)]
    size_list.append(1)


    tree = renumber_tree(var_list)
    trans_list = get_transformer_level1(var_list, tree)

    print("-----orginal tree")
    print_tree(var_list)

    print("-----renumbered tree")
    print_tree(tree)



    max_hash = 1

    for i in range(len(tree)):
        max_hash= max_hash + size_list[i]* max(tree[i])


    return var_index_list,var_len_list, max_hash, size_list, trans_list

def hash_number(y, size_list, ):
    x = 0
    for i in range(len(y)):
        x = x + (y[i] * size_list[i])
    return x


def hash_qname(qname, variational_idx_list,size_list,trans_list ):
    try:

        y = [trans_list[i][qname[variational_idx_list[i]]] for i in range(len(variational_idx_list))]
    except:
        print(len(qname))
    # print(y)
    hash_val = hash_number(y, size_list,)
    return hash_val


def load_hash_chunk(chunk_file):
    df = pd.read_csv(chunk_file)
    chunk_size = df['end'][0] - df['start'][0]
    real_start_list = df['real_start'].to_list()
    hash_rg_list = []
    for i in range(df.shape[0]):
        hash_rg_list.append([ df['start'][i], df['end'][i], df['real_start'][i]])

    return chunk_size,hash_rg_list


def locate_qname(qname,variational_idx_list,size_list,trans_list, chunk_size, hash_rg_list, f_hash_list, pos_width):
    hash_val = hash_qname(qname, variational_idx_list,size_list,trans_list )
    chunk_idx = hash_val // chunk_size
    hash_rg = hash_rg_list[chunk_idx]
    f_hash = f_hash_list[chunk_idx]
    # hash_file = index_dir + "/hashval_%d_%d/merged_gapless.hash"%(hash_rg[0], hash_rg[1])

    real_start = hash_rg[2]

    offset = (hash_val - real_start)*(pos_width+1)

    # f = open(hash_file,'rb')
    f_hash.seek(offset)
    bt = f_hash.read(pos_width)
    # f.close()
    pos = int.from_bytes(bt, 'big')

    return pos


def get_qname_chunk_id(qname,variational_idx_list,size_list,trans_list, chunk_size, hash_rg_list,pos_width ):
    hash_val = hash_qname(qname, variational_idx_list,size_list,trans_list )
    chunk_idx = hash_val // chunk_size 
    hash_rg = hash_rg_list[chunk_idx]
    real_start = hash_rg[2]
    offset = (hash_val - real_start)*(pos_width+1)
    # print(qname, hash_val, chunk_idx, offset)
    return chunk_idx,offset

def extract_read(fqfile,pos_list, qname_list, output_path, qname_separator):
    # pos_list = sort(pos_list)

    with open(output_path,'w') as fw:

        f = open(fqfile,'r')
        for i in tqdm(range(len(pos_list))):
            pos = pos_list[i]
            qname = qname_list[i]
            # print(pos)
            
            f.seek(pos)

            for j in range(4):
                line = f.readline()
                if j == 0:
                    # print(line)
                    qname_detected = line.split()[0].split(qname_separator)[0][1:]
                    # print(qname, qname_detected)
                    assert  qname_detected == qname
                    line = '@' + qname_detected +'\n'
                fw.write(line)
        f.close()

def sort_pos(dc_chunk, dc_chunk_qn):
    for ck_id in dc_chunk:
        pos_list = np.array(dc_chunk[ck_id])
        qn_list = np.array(dc_chunk_qn[ck_id])
        sort_idx = np.argsort(pos_list)
        sorted_pos_list = pos_list[sort_idx]
        sorted_qn_list = qn_list[sort_idx]
        dc_chunk[ck_id] = sorted_pos_list
        dc_chunk_qn[ck_id] = sorted_qn_list
    return dc_chunk, dc_chunk_qn


def locate_one_chunk(ck, hash_rg_list, index_dir, index_pos_list, index_qname_list,
                     pos_width
                     ):
    fq_pos_list= []
    fq_qname_list = []
    # ck = use_ck_list[i]
    hash_rg = hash_rg_list[ck]
    hash_file = index_dir + "/hash/hashval_%d_%d/merged_gapless.hash"%(hash_rg[0], hash_rg[1])
    f_hash = open(hash_file,'rb')
    # index_pos_list = dc_chunk[ck]
    # index_qname_list = dc_chunk_qn[ck]
    for j in tqdm(range(len(index_pos_list)), desc = "process chunk "+str(ck)):
        hash_pos = index_pos_list[j]
        qname = index_qname_list[j]
        f_hash.seek(hash_pos)
        bt = f_hash.read(pos_width)
        pos = int.from_bytes(bt, 'big')
        # print(ck, hash_file, hash_pos, qname, pos)
        fq_pos_list.append(pos)
        fq_qname_list.append(qname)
    f_hash.close()
    # print(fq_pos_list)
    fq_pos_list = np.array(fq_pos_list)
    fq_qname_list = np.array(fq_qname_list)


    #---------sort by fastq pos
    sort_idx = np.argsort(fq_pos_list)
    sorted_fq_pos_list = fq_pos_list[sort_idx]
    sorted_fq_qname_list = fq_qname_list[sort_idx]
    return sorted_fq_pos_list, sorted_fq_qname_list

def cut_2_lists(list1,list2,n_ck):
    size = int(len(list1)/n_ck) + 1
    new_list1 = []
    new_list2 = []
    for i in range(0, len(list1), size):
        list1_seg = list1[i:i+size]
        list2_seg = list2[i:i+size]
        new_list1.append(list1_seg)
        new_list2.append(list2_seg)
    return new_list1, new_list2

def merge_files(inlist, outfile):
    cmd = "cat " + " ".join(inlist)+ " > " + outfile
    Popen(cmd, shell = True).wait()

    cmd = "rm " + " ".join(inlist)
    Popen(cmd, shell = True).wait()


def retrieve_reads_multi_threads(fastq_file, qname_file, output_path, index_dir, n_thread, qname_separator):
    stats_pickle = index_dir+"/read_name_stats.p"
    with open(stats_pickle,'rb') as f:
        dc = pickle.load(f)
    print("read name lengths:", list(dc.keys()))
    if len(dc) >1:
        print("Can only handle equal length read names so far")
        exit()

    qname_len = list(dc.keys())[0]
    var_path = index_dir+f"/read_name_stats_len{qname_len}_all.txt"

    var_index_list,var_len_list, max_hash, size_list, trans_list = load_var(var_path)
    chunk_size,hash_rg_list = load_hash_chunk(index_dir + "/hash/chunks_hash_realst.csv")

    hash_stats_path = index_dir+"/hash/hash_stats.json"
    with open(hash_stats_path,'r') as f:
        hash_dc = eval(f.read())
    pos_width = hash_dc['pos_width']

    # ---- hash and locate all qnames
    with open(qname_file,'r') as f:
        qname_list = f.read().split('\n')[:-1]
    dc_chunk = defaultdict(list)
    dc_chunk_qn = defaultdict(list)
    for qname in tqdm(qname_list, desc = "hash qnames"):
        ck_id,pos = get_qname_chunk_id(qname,var_index_list,size_list,trans_list, chunk_size, hash_rg_list,pos_width )
        dc_chunk[ck_id].append(pos)
        dc_chunk_qn[ck_id].append(qname)

    #---------sort by index pos
    dc_chunk, dc_chunk_qn = sort_pos(dc_chunk, dc_chunk_qn)
    # print(dc_chunk)

    use_ck_list = sorted(list(dc_chunk.keys()))
    print("num chunks:",len(dc_chunk))

    for i in range(len(use_ck_list)):
        ck = use_ck_list[i]
        # print(dc_chunk[ck])
        print(f"{i+1} chunk ",ck,': ', len(dc_chunk[ck]),' reads')

    # ---- locate all qnames in fastq


    sequences = Parallel(n_jobs=n_thread)(delayed(locate_one_chunk )\
                                        (ck, hash_rg_list, index_dir, \
                                        dc_chunk[ck], dc_chunk_qn[ck],
                        pos_width) for ck in use_ck_list)

    fq_pos_list, fq_qname_list = [], []

    for sorted_fq_poss, sorted_fq_qnames in sequences:
        fq_pos_list.extend(sorted_fq_poss)
        fq_qname_list.extend(sorted_fq_qnames)

    #---------sort by fastq pos
    fq_pos_list = np.array(fq_pos_list)
    fq_qname_list = np.array(fq_qname_list)
    sort_idx = np.argsort(fq_pos_list)
    sorted_fq_pos_list = fq_pos_list[sort_idx]
    sorted_fq_qname_list = fq_qname_list[sort_idx]

    extract_read(fastq_file,sorted_fq_pos_list, sorted_fq_qname_list, output_path, qname_separator)


    # fq_pos_seg_list, fq_qname_seq_list = cut_2_lists(sorted_fq_pos_list,sorted_fq_qname_list,n_thread)
    # out_list = [output_path+"."+str(i+1) for i in range(len(fq_pos_seg_list))]



    # sequences = Parallel(n_jobs=n_thread)(delayed(extract_read)\
    #                                       (fastq_file,fq_pos_seg_list[i], fq_qname_seq_list[i], out_list[i])\
    #                                       for i in range(len(fq_pos_seg_list)))
    

    # merge_files(out_list, output_path)
                                          




if __name__ == "__main__":
    retrieve_reads_multi_threads(fastq_file, qname_file, output_path, index_dir, n_thread, qname_separator)