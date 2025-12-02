import argparse
from argparse import ArgumentParser
parser = ArgumentParser(description="",
	usage='use "python3 %(prog)s --help" for more information',
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--indir','-i')
parser.add_argument('--chr_fastq','-fq')
# parser.add_argument('--read_type','-rt', choices = ['PE','SE'])
parser.add_argument('--max_mem','-mem', help = "GB", default= 100, type = int)
parser.add_argument('--chr_num','-chr', type = int )
parser.add_argument('--n_thread','-t', type = int, default = 22 )
parser.add_argument('--sample_name','-name',help="Required parameter; sample name you can define, for example, S12878",required=True)
parser.add_argument('--delete_temp_file','-d', action='store_true')
args = parser.parse_args()
indir = args.indir
# read_type = args.read_type
max_mem = args.max_mem
# output_path = args.output_path
n_thread = args.n_thread
chr_fastq = args.chr_fastq
chr_num = args.chr_num
sample_name = args.sample_name

import pdb
#pdb.set_trace()
from collections import defaultdict
import pickle
import os
from argparse import ArgumentParser
from tqdm import tqdm
import glob
from subprocess import Popen
import psutil
from datetime import datetime
from joblib import Parallel, delayed
import pandas as pd
import psutil, os, gc

def save_pickle_file(dict1,fname):
    for value in dict1:
        dict1[value] = dict(dict1[value])
    my_dict = dict(dict1)
    with open(fname, "wb") as f:
        pickle.dump(my_dict,f) 


def Extract_mole_num_from_phased_mole(phased_file,PS_flag_dict_cut_file,chr_num):
    phase_dict = defaultdict(list)
    PS_flag_dict_cut = pickle.load(open(PS_flag_dict_cut_file,"rb"))
    f = open(phased_file,"r")
    for line in f:
        data = line.rsplit()
        PS_flag_info = data[-2]
        PS_flag = int(PS_flag_info.split(":")[1])
        HP_flag = data[-1]
        mole_num = int(data[6])
        PS_flag_end = PS_flag_dict_cut[PS_flag]
        phase_dict[(PS_flag,PS_flag_end,HP_flag)].append(mole_num)

    #pickle.dump(phase_dict,open("phase_dict_chr" + str(chr_num) + ".p","wb"))
    print("done~")
    return phase_dict

def check_qname_in_PS(block_info,qname_pos_list):
    block_start = block_info[0]
    block_end = block_info[1]
    flag_in = 0
    for pos in qname_pos_list:
        if pos <= block_end and pos >= block_start:
            flag_in = 1
            break
    return flag_in 


def Extract_qname(phased_dict_mole_num,mole_qname_dict,qname_pos_dict,barcoded_fastq_file,chr_num,output_dir,):
    # assert read_type in ['SE','PE']

    phased_dict_qname_2 = {}
    curr = 0
    flag_in = 0
    for key, mole_num_list in tqdm(phased_dict_mole_num.items(), desc = " extract qname - hp dict"):
        curr += 1
        for mole_num in mole_num_list:
            qname_list = mole_qname_dict[mole_num]
            for qname in qname_list:
                flag_in = check_qname_in_PS(key,qname_pos_dict[qname])
                if flag_in == 1:
                    phased_dict_qname_2[qname] = key
    if not os.path.exists(output_dir):
        os.system("mkdir -p " + output_dir)
    outfile = output_dir + "/phased_dict_qname_2.p"
    # save_pickle_file(phased_dict_qname_2, outfile)
    with open(outfile,'wb') as fw:
        pickle.dump(phased_dict_qname_2, fw)
    

def split_fastq_by_lines(fastq_file, out_dir):
    split_cmd = "cat " + fastq_file + " | split -l 5000000 - " + out_dir + "/barcoded.fastq_part"
    print(split_cmd)
    Popen(split_cmd,shell=True).wait()
    fq_list = glob.glob(out_dir + "/barcoded.fastq_part*")
    return fq_list



def split_one_fq(barcoded_fastq_file,phased_dict_qname_2_path, out_dir, start_byte, end_byte):
    phased_dict_qname_2 = pickle.load(open(phased_dict_qname_2_path,'rb'))
    '''
    +
    qual
    @
    seq
    +\n

    0 or 4 : 2

    qual
    @
    seq
    +\n
    qual
    3:1



    @
    seq
    +\n
    qual
    @
    2:0 or 4


    seq
    +\n
    qual
    @
    seq
    1:3
    '''

    f = open(barcoded_fastq_file,"r")
    if start_byte!=0:
        f.seek(start_byte)
        # f.readline()
        found_plus = 0
        byte_list = []
        lines = []
        for i in range(5):
            byte_position = f.tell()   # record current byte offset
            byte_list.append(byte_position)
            line = f.readline()
            lines.append(line)
            if line == '+\n':
                plus_i = i
                found_plus=1
                # break 
        if found_plus == 0:
            print(f"could not find '+\\n' line after searching 5 rows starting from {start_byte}, please check your fastq file")
            exit()
        if plus_i in [0,4]:
            start_i = 2
        elif plus_i == 3:
            start_i = 1
        elif plus_i == 2:
            if lines[0][0]=='@':
                start_i = 0 
            else:
                start_i = 4 
        else:
            start_i = 3
        
        actual_start_byte = byte_list[start_i]
    else:
        actual_start_byte = start_byte

    # actual_start_byte = f.tell()
    f.seek(actual_start_byte)
    # then it will start with read name line
    dc = defaultdict(list)
    count = 0
    # f.seek(real_start_byte) 
    while True:
        byte_position = f.tell()   # record current byte offset
        line = f.readline()
        count+=1
        if not line:
            break
        data = line.rsplit()  # decode if needed
        qname_curr = data[0][1:]
        if qname_curr in phased_dict_qname_2:
            _PS_HP = phased_dict_qname_2[qname_curr]
            dc[_PS_HP].append(byte_position)
        for i in range(3):
            f.readline()
            count += 1
        if (count!=1) & (count % 1000000 == 0):
            print(f"processed {barcoded_fastq_file} {count} lines...")
        cur_byte = f.tell()
        if cur_byte >= end_byte:
            break
    actual_end_byte = cur_byte
    outfile = out_dir + f'/HP_Reads_{start_byte}_{end_byte}.p'
    with open(outfile, "wb") as f:
        pickle.dump(dc,f)

    outlog = out_dir + f'/HP_Reads_{start_byte}_{end_byte}.log'
    with open(outlog, "w") as f:
        f.write(f"original start byte {start_byte}, original end byte {end_byte}\n")
        f.write(f"actual start byte {actual_start_byte}, actual end byte {actual_end_byte}\n")

    print("finished extracting " + barcoded_fastq_file + f'({start_byte}-{end_byte}), saved to '+ outfile)
    return outfile, actual_start_byte, actual_end_byte


def estimate_memory_threads(file_path, max_mem, log_path):
    """
    Load a pickled dictionary, compute memory usage delta (MB), and log to file.
    """
    process = psutil.Process(os.getpid())
    mem_before = process.memory_info().rss

    # Load the pickled dict
    with open(file_path, "rb") as f:
        data = pickle.load(f)

    mem_after = process.memory_info().rss
    mem_used_gb = (mem_after - mem_before) / (1024 ** 3)

    

    max_t = int(max_mem/(mem_used_gb*1.5))

    # Write to log file with timestamp
    with open(log_path, "w") as log:
        log.write(f"{datetime.now():%Y-%m-%d %H:%M:%S} | {file_path} | {mem_used_gb:.2f} GB\n")
        log.write(f"max allowed mem: {max_mem}GB\n")
        log.write(f"max allowed threads: {max_t}\n")

    print(f"Loaded {file_path}")
    print(f"Memory usage increase: {mem_used_gb:.2f} GB")
    print(f"max allowed threads: {max_t}\n")
    print(f"Log saved to {log_path}")

    return mem_used_gb, max_t

def split_file_into_chunks(file, n_chunk, outfile):
    size = os.path.getsize(file)
    print(size)
    chunk_size = size// n_chunk
    
    range_list = [ (i*chunk_size, (i+1) * chunk_size) for i in range(n_chunk-1)]
    range_list.append((( n_chunk-1)* chunk_size, size))
    df = pd.DataFrame(range_list, columns=['start_byte', 'end_byte'])
    df.to_csv(outfile, sep = '\t')
    return range_list

def validate(chunks, sequences, log_file):
    cnt = 0
    f = open(log_file,'w')
    for i, chunk in enumerate(chunks[:-1]):
        cur_end  = sequences[i][2]
        next_start  = sequences[i+1][1]
        if cur_end!=next_start:
            msg = f"chunk {str(chunks[i])} and {str(chunks[i+1])} not consecutive! {cur_end}!={next_start}"
            print(msg)
            f.write(msg)
            cnt+=1
    if cnt ==0:
        msg = "All consecutive!"
        print(msg)
        f.write(msg)
        f.close()
    else:
        print("Did not pass validation! Check your block!")
        f.close()
        exit()

def merge_result(out_dir, chunks, out_file):
    files = [out_dir + f'/HP_Reads_{chunk[0]}_{chunk[1]}.p'  for chunk in chunks]
    dc = pickle.load(open(files[0],'rb'))
    for file in tqdm(files[1:], desc = "merging results"):
        dc_cur = pickle.load(open(file,'rb'))
        for key, val in dc_cur.items():
            if key in dc:
                prev_val = dc[key]
                try:
                    assert prev_val[-1] < val[0]
                except:
                    print(f"found error in file {file} with previous blocks; PS_tag {key} previous end {prev_val[-1]} >= current start {val[0]}")
                    exit()
                dc[key].extend(val)
            else:
                dc[key] = val 
    with open(out_file,'wb') as f:
        pickle.dump(dc, f)

def split_all_fq( fq_file,  output_dir, n_thread, max_mem):

    phased_dict_qname_2_path = output_dir + "/phased_dict_qname_2.p"
    log_path = phased_dict_qname_2_path+"_mem_estimation.log"

    dc_mem, max_t = estimate_memory_threads(phased_dict_qname_2_path, max_mem, log_path)
    # max_t = 15
    use_thread = min(max_t, n_thread)
    with open(log_path,'a') as f:
        f.write(f"actually used threads: {use_thread}\n")
    chunks = split_file_into_chunks(fq_file, 30,  output_dir + "/fastq_split.tsv")
    sequences = Parallel(n_jobs=use_thread)(delayed(split_one_fq)(fq_file, phased_dict_qname_2_path,  output_dir, chunk[0], chunk[1]) for chunk in tqdm(chunks, desc = "process by block"))
    validate(chunks, sequences, output_dir+"/validate.log")
    merge_result(output_dir, chunks, output_dir+"/HP_Reads_merged.p")






def clean_folder(output_dir):
    cmd = f"rm {output_dir}/HP_Reads_*_*.* {output_dir}/phased_dict_qname_2.p"
    print(cmd)
    Popen(cmd, shell=True).wait()








def Extract_qname_hp_dict(phased_h5_file,PS_flag_dict_cut_file,chr_num,mole_qname_dict_file,output_dir,logger):

    process = psutil.Process(os.getpid())

    mem0=process.memory_info().rss / 1024**2
    phased_dict_mole_num = Extract_mole_num_from_phased_mole(phased_h5_file,PS_flag_dict_cut_file,chr_num)
    mem1=process.memory_info().rss / 1024**2
    logger.info(f"Estimated phased_dict_mole_num memory usage: {(mem1-mem0):.2f} MB")

    mem0 = process.memory_info().rss / 1024**2
    mole_qname_dict = pickle.load(open(mole_qname_dict_file,"rb"))
    mem1 = process.memory_info().rss / 1024**2
    logger.info(f"Estimated mole_qname_dict memory usage: {(mem1-mem0):.2f} MB")

    phased_dict_qname_2 = {}
    curr = 0
    flag_in = 0
    for key, mole_num_list in tqdm(phased_dict_mole_num.items(), desc = " extract qname - hp dict"):
        curr += 1
        for mole_num in mole_num_list:
            qname_dict = mole_qname_dict[mole_num]
            for qname in qname_dict:
                flag_in = check_qname_in_PS(key,qname_dict[qname])
                if flag_in == 1:
                    phased_dict_qname_2[qname] = key
    if not os.path.exists(output_dir):
        os.system("mkdir -p " + output_dir)
    outfile = output_dir + "/phased_dict_qname_2.p"
    # save_pickle_file(phased_dict_qname_2, outfile)
    with open(outfile,'wb') as fw:
        pickle.dump(phased_dict_qname_2, fw)

def Extract_start(output_dir,chr_num,phased_h5_file,PS_flag_dict_cut_file,mole_qname_dict_file,qname_pos_dict_file,chr_fastq,max_mem, n_thread):
    import logging
    ## set logger
    logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')
    logger = logging.getLogger(" ")

    process = psutil.Process(os.getpid())
    # logger.info("0")
    # logger.info(f"current memory usage: {process.memory_info().rss / 1024**2:.2f} MB")
    # phased_dict_mole_num = Extract_mole_num_from_phased_mole(phased_h5_file,PS_flag_dict_cut_file,chr_num)
    # logger.info("1")
    # logger.info(f"current memory usage: {process.memory_info().rss / 1024**2:.2f} MB")
    # mole_qname_dict = pickle.load(open(mole_qname_dict_file,"rb"))
    # logger.info("2")
    # logger.info(f"current memory usage: {process.memory_info().rss / 1024**2:.2f} MB")
    # qname_pos_dict = pickle.load(open(qname_pos_dict_file,"rb"))
    # logger.info("3")
    # logger.info(f"current memory usage: {process.memory_info().rss / 1024**2:.2f} MB")
    # Extract_qname(phased_dict_mole_num,mole_qname_dict,qname_pos_dict,chr_fastq,chr_num,output_dir,)
    # logger.info("4")
    # logger.info(f"current memory usage: {process.memory_info().rss / 1024**2:.2f} MB")
    

    
    # logger.info(f"Before cleanup: {process.memory_info().rss / 1024**2:.2f} MB")

    # del phased_dict_mole_num, mole_qname_dict, qname_pos_dict
    # gc.collect()

    # logger.info(f"After cleanup: {process.memory_info().rss / 1024**2:.2f} MB")

    Extract_qname_hp_dict(phased_h5_file,PS_flag_dict_cut_file,chr_num,mole_qname_dict_file,output_dir,logger)
    split_all_fq( chr_fastq,  output_dir, n_thread, max_mem)
    clean_folder(output_dir)






# chr_num=22
phased_h5_file=indir + f"/phase_blocks_cut_highconf/chr{chr_num}.phased_final_cut_by_100000"    
PS_flag_dict_cut_file= indir + f"/phase_blocks_cut_highconf/chr{chr_num}.phased_final_cut_by_100000_phase_blocks.p"
mole_qname_dict_file= indir +  f"/H5_for_molecules/{sample_name}_chr{chr_num}_qname.p"
qname_pos_dict_file= indir + f"/H5_for_molecules/{sample_name}_chr{chr_num}_qname_pos.p"
output_dir = indir + f'/Local_Assembly_by_chunks/chr{chr_num}_files_cutPBHC/'
# chr_fastq="/data/maiziezhou_lab/Datasets/stLFR_data/NA24385_giab/fastq_by_chr/chr22/NA24385_stlfr_giab_chr22.fastq"
# read_type = "PE"
# n_thread = 10


Extract_start(output_dir,chr_num,phased_h5_file,PS_flag_dict_cut_file,
              mole_qname_dict_file,qname_pos_dict_file,chr_fastq,max_mem, n_thread)
