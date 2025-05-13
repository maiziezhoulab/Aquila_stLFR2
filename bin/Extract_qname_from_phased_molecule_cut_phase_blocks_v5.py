import argparse
from argparse import ArgumentParser
parser = ArgumentParser(description="",
	usage='use "python3 %(prog)s --help" for more information',
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--indir','-i')
parser.add_argument('--chr_fastq','-fq')
parser.add_argument('--read_type','-rt', choices = ['PE','SE'])
parser.add_argument('--chr_num','-chr', type = int )
parser.add_argument('--n_thread','-t', type = int, default = 22 )
parser.add_argument('--sample_name','-name',help="Required parameter; sample name you can define, for example, S12878",required=True)
parser.add_argument('--delete_temp_file','-d', action='store_true')
args = parser.parse_args()
indir = args.indir
read_type = args.read_type
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

from joblib import Parallel, delayed



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


def Extract_qname(phased_dict_mole_num,mole_qname_dict,qname_pos_dict,barcoded_fastq_file,chr_num,output_dir,read_type):
    assert read_type in ['SE','PE']

    phased_dict_qname_2 = defaultdict(lambda: defaultdict(int))
    qname_total_dict = defaultdict(int)
    curr = 0
    flag_in = 0
    for key, mole_num_list in tqdm(phased_dict_mole_num.items(), desc = " extract qname - hp dict"):
        curr += 1
        for mole_num in mole_num_list:
            qname_list = mole_qname_dict[mole_num]
            for qname in qname_list:
                flag_in = check_qname_in_PS(key,qname_pos_dict[qname])
                if flag_in == 1:
                    phased_dict_qname_2[qname][key] = 1
    if not os.path.exists(output_dir):
        os.system("mkdir -p " + output_dir)
    outfile = output_dir + "/phased_dict_qname_2.p"
    save_pickle_file(phased_dict_qname_2, outfile)
    

def split_fastq_by_lines(fastq_file, out_dir):
    split_cmd = "cat " + fastq_file + " | split -l 5000000 - " + out_dir + "/barcoded.fastq_part"
    print(split_cmd)
    Popen(split_cmd,shell=True).wait()
    fq_list = glob.glob(out_dir + "/barcoded.fastq_part*")
    return fq_list



def split_one_fq(barcoded_fastq_file,phased_dict_qname_2_path, read_type, out_dir):
    phased_dict_qname_2 = pickle.load(open(phased_dict_qname_2_path,'rb'))

    


    f = open(barcoded_fastq_file,"r")
    count = 0
    flag = 0

    if read_type == 'SE':
        unit = 4

    else:
        unit = 8


    dc = defaultdict(list)

    for line in f:
        
        data = line.rsplit()
        if count%unit == 0:
            qname_curr = data[0][1:]
            if qname_curr in phased_dict_qname_2:
                flag = 1
                PS_flag_info = phased_dict_qname_2[qname_curr]
                for _PS_HP, _val in PS_flag_info.items():
                    _PS_flag = _PS_HP[0]
                    _PS_flag_end = _PS_HP[1]
                    _HP_flag = _PS_HP[2]
                file_curr = output_dir +"fastq_by_" + str(_PS_flag) + "_" + str(_PS_flag_end) + "_" + _HP_flag + ".fastq"
                dc[file_curr].append(line)
                # if os.path.isfile(file_curr):
                #     fw_curr = open(file_curr,"a")
                #     fw_curr.write(line)
                # else:
                #     fw_curr = open(file_curr,"w")
                #     fw_curr.write(line)
                del phased_dict_qname_2[qname_curr]

            # else:
            #     del phased_dict_qname_2[qname_curr]

        elif count%unit == (unit-1):
            if flag == 1:
                flag = 0

                dc[file_curr].append(line)

                # fw_curr.write(line)
                # fw_curr.close()

        else:
            if flag == 1:
                dc[file_curr].append(line)
                # fw_curr.write(line)
        


        count += 1
        if (count!=1) & (count % 1000000 == 0):
            print(f"processed {barcoded_fastq_file} {count} lines...")

    outfile = out_dir + '/HP_Reads_'+barcoded_fastq_file.split('/')[-1]+'.p'

    with open(outfile, "wb") as f:
        pickle.dump(dc,f) 


    print("finished extracting " + barcoded_fastq_file + ', saved to '+ outfile)
    return outfile


def split_all_fq( read_type, output_dir, n_thread):
    fq_list = glob.glob(output_dir + "/barcoded.fastq_part*")
    phased_dict_qname_2_path = output_dir + "/phased_dict_qname_2.p"

    sequences = Parallel(n_jobs=n_thread)(delayed(split_one_fq)(fq, phased_dict_qname_2_path, read_type, output_dir) for fq in tqdm(fq_list))


def merge_fq(output_dir):
    dc_file_list = glob.glob(output_dir + "/HP_Reads_barcoded.fastq_part*.p")
    dc_final = defaultdict(list)

    for dc_file in tqdm(dc_file_list, desc = "reduce result"):
        dc_new = pickle.load(open(dc_file,'rb'))

        for key, val in dc_new.items():
            dc_final[key].extend(val)

        del dc_new 

    for key, val in tqdm(dc_final.items(), desc = "write fastqs"):
        with open(key, 'w') as fw:
            fw.writelines(val)

        # del dc_final[key]



def clean_folder(output_dir):
    cmd = f"rm {output_dir}/HP_Reads_barcoded.fastq_part*.p {output_dir}/barcoded.fastq_part*"
    print(cmd)
    Popen(cmd, shell=True).wait()










def Extract_start(output_dir,chr_num,phased_h5_file,PS_flag_dict_cut_file,mole_qname_dict_file,qname_pos_dict_file,chr_fastq,read_type, n_thread):
    phased_dict_mole_num = Extract_mole_num_from_phased_mole(phased_h5_file,PS_flag_dict_cut_file,chr_num)
    mole_qname_dict = pickle.load(open(mole_qname_dict_file,"rb"))
    qname_pos_dict = pickle.load(open(qname_pos_dict_file,"rb"))
    
    Extract_qname(phased_dict_mole_num,mole_qname_dict,qname_pos_dict,chr_fastq,chr_num,output_dir,read_type)
    split_fastq_by_lines(chr_fastq, output_dir)
    split_all_fq(  read_type, output_dir, n_thread)
    merge_fq(output_dir)
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
              mole_qname_dict_file,qname_pos_dict_file,chr_fastq,read_type, n_thread)
