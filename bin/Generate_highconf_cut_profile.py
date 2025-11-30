#!/usr/bin/env python
import os
import numpy as np
from argparse import ArgumentParser
from multiprocessing import Pool,cpu_count
import pickle
from collections import defaultdict
import subprocess
import pandas as pd 
from datetime import datetime
from tqdm import tqdm
import gc
import psutil


def get_read_depth(bam_file,mbq,mmq,chr_num,out_dir,threads, logger):
    if chr_num == 23:
        chr_num = "X"
    else:
        chr_num = str(int(chr_num))
    if not os.path.exists(out_dir):
        os.system("mkdir -p " + out_dir)
    depth_file_noHclip = out_dir + "/chr" + chr_num + "_depth_noHclip.txt"

    depth_cmd = "samtools view -b -u -h -e 'cigar!~\"[HS]\"' -@ "+ str(threads) +" "+ bam_file + " chr" + chr_num + " | " \
        + "samtools depth -q " + str(mbq) +  " -Q " + str(mmq) + " - > " + depth_file_noHclip

    logger.info(depth_cmd)
    subprocess.Popen(depth_cmd,shell=True).wait()


    return depth_file_noHclip 
    






def get_global_track_for_breakpoints(depth_file_noHclip , uniq_map_file, output_file, logger):
    process = psutil.Process(os.getpid()) 
    logger.info(f"step0.1 current  Memory usage: {process.memory_info().rss / 1024 ** 2:.2f} MB")
    logger.info(f"Step0 starts. {datetime.now()}")
    depth_info_temp = pd.read_csv(depth_file_noHclip, sep = '\t', header = None).iloc[:,1:].values
    logger.info(f"step0.2 current  Memory usage: {process.memory_info().rss / 1024 ** 2:.2f} MB")
    depth_info_temp[:,0] = depth_info_temp[:,0]+1
    logger.info(f"Step0 finished. {datetime.now()}")
    logger.info(f"step0.3 current  Memory usage: {process.memory_info().rss / 1024 ** 2:.2f} MB")


    depth_info_temp = np.array(depth_info_temp, dtype=int)
    logger.info(f"step0.4 current  Memory usage: {process.memory_info().rss / 1024 ** 2:.2f} MB")

    # Convert to NumPy array once
    


    # Compute mean and threshold
    mean_cov = depth_info_temp[:, 1].mean()
    logger.info(f"step0.5 current  Memory usage: {process.memory_info().rss / 1024 ** 2:.2f} MB")
    cov_threshold = mean_cov * 0.8
    logger.info(f"mean coverage: {mean_cov:.2f}, coverage thresh: {cov_threshold:.2f}")

    # --- Step 1: filter coverage first (cheap vectorized filter)
    depth_info = depth_info_temp[depth_info_temp[:, 1] > cov_threshold,:]
    logger.info(f"Step1 finished. {datetime.now()}")
    logger.info(f"step1.1 current  Memory usage: {process.memory_info().rss / 1024 ** 2:.2f} MB")
    del depth_info_temp
    gc.collect()
    np.empty(0)  # forces some allocator cleanup
    logger.info(f"step1.2 current  Memory usage (supposed to drop): {process.memory_info().rss / 1024 ** 2:.2f} MB")

    # --- Step 2: intersect using NumPy for speed ---
    pos_set = np.array(list(set(depth_info[:, 0].flatten())))
    logger.info(f"step2.0 current  Memory usage: {process.memory_info().rss / 1024 ** 2:.2f} MB")
    del depth_info
    gc.collect()
    logger.info(f"step2.1 current  Memory usage (supposed to drop?): {process.memory_info().rss / 1024 ** 2:.2f} MB")
    with open(uniq_map_file, "rb") as f:
        uniq_map = np.array(list(pickle.load(f).keys()))
        logger.info(f"step3.1 current  Memory usage: {process.memory_info().rss / 1024 ** 2:.2f} MB")



    pos_set.sort(); uniq_map.sort()
    logger.info(f"step3.2 current  Memory usage: {process.memory_info().rss / 1024 ** 2:.2f} MB")
    pos_set = np.intersect1d(pos_set, uniq_map , assume_unique=True)

    logger.info(f"step4.1 current  Memory usage: {process.memory_info().rss / 1024 ** 2:.2f} MB")
    del uniq_map
    gc.collect()
    logger.info(f"step4.2 current  Memory usage (supposed to drop): {process.memory_info().rss / 1024 ** 2:.2f} MB")
    logger.info(f"# valid postions:{len(pos_set)} ")
    logger.info(f"Step2 finished. {datetime.now()}")

    # arr = np.array(list(pos_set),  dtype=np.int32)
    np.save(output_file, pos_set)




def gen_highconf_profile(bam_file,chr_start,chr_end,mbq,mmq,out_dir,uniq_map_dir,threads):
    import logging
    ## set logger
    logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')
    logger = logging.getLogger(" ")


    for chr_num in range(chr_start,chr_end + 1):
        print("Processing chr" + str(chr_num))
        process = psutil.Process(os.getpid())
        before_mem =  process.memory_info().rss / 1024 ** 2
        depth_info_file = get_read_depth(bam_file, mbq, mmq, chr_num, out_dir, threads, logger)    
        after_mem =  process.memory_info().rss / 1024 ** 2
        logger.info(f"[depth info]  Memory usage: {after_mem-before_mem:.2f} MB")
        logger.info(f"current  Memory usage: {after_mem:.2f} MB")
        uniq_map_file = uniq_map_dir + "/uniq_map_chr" + str(chr_num) + ".p"
        
        output_file = out_dir + "chr" + str(chr_num) + "_global_track.npy"
        get_global_track_for_breakpoints(depth_info_file,uniq_map_file,output_file, logger)
    print("Finishing generating highconf profile")

    cmd = f"rm {out_dir}/chr*_depth_noHclip.txt"
    os.system(cmd)


if __name__ == "__main__":
    parser = ArgumentParser(description="Run depth all:")
    parser.add_argument('--bam_file','-bam',help="input sorted bam file")
    parser.add_argument('--chr_start','-start',type=int,help="chromosome start from", default=1)
    parser.add_argument('--chr_end','-end',type=int,help="chromosome end by", default=23)
    parser.add_argument('--min_base_qual','-minbq',type=int,help="minimum base quality to consider a base for haplotype fragment, default 13", default=13)
    parser.add_argument('--min_map_qual','-minmq',type=int,help="minimum read mapping quality to consider a read for phasing, default 20", default=20)
    parser.add_argument('--out_dir','-o_dir', help="Directory to store outputs")
    parser.add_argument('--uniq_map_dir','-uniq_dir', help="Directory to 100-mer uniqness")
    parser.add_argument('--threads','-t',type=int,help="number of threads for running, default 20", default=20)
    args = parser.parse_args()

    out_dir = args.out_dir
    if os.path.exists(out_dir):
        print("using existing output folder: " + out_dir)
    else:
        os.makedirs(out_dir)

    bam_file = args.bam_file
    chr_start = args.chr_start
    chr_end = args.chr_end
    mbq = args.min_base_qual
    mmq = args.min_map_qual
    uniq_map_dir = args.uniq_map_dir
    num_threads = args.threads

    gen_highconf_profile(bam_file,chr_start,chr_end,mbq,mmq,out_dir,uniq_map_dir,num_threads)
