





from subprocess import Popen
import os
from load_read1 import retrieve_reads
from joblib import Parallel, delayed
from tqdm import tqdm

def load_names(nf):
    with open(nf,'r') as f:
        names  = f.read().split('\n')[:-1]
    return names

def cut_chunks(names,n_ck, outdir, prefix):
    size_per_ck = int(len(names) / n_ck) + 1
    chunks  = [(i*size_per_ck, (i+1)*size_per_ck) for i in range(n_ck-1)]
    chunks.append((size_per_ck*(n_ck-1), len(names)))

    outfile = outdir+'/chunks_'+prefix+'.txt'
    with open(outfile,'w') as f:
        print(chunks, file = f)
    return chunks 

def merge(in_list, out_file):
    cmd = 'cat '+ ' '.join(in_list) + " > " + out_file
    print(cmd)
    Popen(cmd, shell= True).wait()
    return 


def extract_fq_from_fq(out_dir,n_ck,fqfile, index_dir, prefix,n_thread):
    # n_ck = 70
    unused_name_uniq  = out_dir+'/unused_name_uniq_sorted.txt'
    unsed_names  = load_names(unused_name_uniq)
    chunks = cut_chunks(unsed_names, n_ck, out_dir, prefix)
    print(len(unsed_names))
    print(chunks)

    outfile_list = [  out_dir+f'/{prefix}_{ck[0]}_{ck[1]}.fq' for ck in chunks]
    sequences = Parallel(n_jobs=n_thread)(delayed(retrieve_reads)(unsed_names[chunks[i][0]: chunks[i][1]],index_dir,fqfile, outfile_list[i])  for i in tqdm(range(len(chunks))))
    merged_fq = out_dir+f'/merged_{prefix}.fq'
    merge(outfile_list, merged_fq)
    return 


if __name__ == '__main__':
    import argparse
    from argparse import ArgumentParser
    parser = ArgumentParser(description="",
        usage='use "python3 %(prog)s --help" for more information',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--fq_file','-fq')
    parser.add_argument('--index_dir','-index')
    parser.add_argument('--out_dir','-o')
    parser.add_argument('--n_chunk','-nck', type = int, help = "n chunks")

    parser.add_argument('--prefix','-px',)
    parser.add_argument('--n_thread','-t', type = int, default = 22 )

    args = parser.parse_args()
    index_dir = args.index_dir
    fqfile = args.fq_file
    out_dir = args.out_dir
    n_ck = args.n_chunk
    prefix = args.prefix
    n_thread = args.n_thread
    extract_fq_from_fq(out_dir,n_ck,fqfile, index_dir, prefix,n_thread)