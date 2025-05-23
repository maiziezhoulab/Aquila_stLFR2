import pdb
#pdb.set_trace()
import pickle
import gzip
import glob
from collections import defaultdict
import numpy as np
from subprocess import Popen
from multiprocessing import Pool,cpu_count,active_children
import time
import os
import sys
from argparse import ArgumentParser
parser = ArgumentParser(description="Author: xzhou15@cs.stanford.edu\n",usage='use "python3 %(prog)s --help" for more information')
parser.add_argument('--assembly_vcf','-v_assembled',help="vcf file",required=True)
parser.add_argument('--vcf_file','-v',help="vcf file",required=True)
parser.add_argument('--chr_start','-start',type=int,help="chromosome start from", default=1)
parser.add_argument('--chr_end','-end',type=int,help="chromosome end by", default=23)
parser.add_argument('--out_dir','-o', help="Directory to store output",default="./results_phasing_vcf")
parser.add_argument('--assembly_dir','-a_dir', help="folder to store assembled results",required=True)
parser.add_argument('--block_len_use','-bl',type=int,help="phase block len threshold",default=100000)
args = parser.parse_args()
__author__ = "Xin Zhou@Stanford"


# 1. Phasing assembled variants by contigs and phase blocks
def check_PS_flag_in_HCbk(HCbk_dict,PS_flag_1,PS_flag_2):
    flag = 0
    flag_1 = 0
    flag_2 = 0
    for key,val in HCbk_dict.items():
        flag_1 = int(key[0])
        flag_2 = int(key[1])
        if PS_flag_1 >= flag_1 and PS_flag_2 <= flag_2:
            flag = 1
            break

    return (flag,flag_1,flag_2)


def Phasing_assembled_SV(HCbk_dict,vcf_file,chr_num,out_dir):
    phased_variants = defaultdict(lambda: defaultdict(list))
    phased_variants_loci = defaultdict(lambda: defaultdict(list))
    #f = gzip.open(vcf_file,"rb")
    f = open(vcf_file,"rb")
    if chr_num == 23:
        chr_num = "X"
    for line in f:
        data = line.decode().rsplit()
        if data[0][0] != "#":
            chr_num_curr = data[0]
            if chr_num_curr == "chr" + str(chr_num):
                GT = data[9].split(":")[0]
                if GT == "0/1":
                    if int(data[9].split(":")[1].split("_")[0]) == 0:
                        PS_info = data[9].split("PS")[2].split("_")
                        PS_flag_1 = int(PS_info[0])
                        PS_flag_2 = int(PS_info[1])
                        hp_flag = PS_info[2]
                    else:
                        PS_info = data[9].split("PS")[1].split("_")
                        PS_flag_1 = int(PS_info[0])
                        PS_flag_2 = int(PS_info[1])
                        hp_flag = PS_info[2]
                    
                    ## check if PS_flag_1 is in certain region or not?
                    flag,flag_1,flag_2 = check_PS_flag_in_HCbk(HCbk_dict,PS_flag_1,PS_flag_2)

                    if flag == 1:
                        # use the start and end in Hcbk dict key
                        phased_variants[(flag_1,flag_2)][hp_flag].append(data)
                        phased_variants_loci[(flag_1,flag_2)][hp_flag].append(int(data[1]))
                    else:
                        # use the original PS tag
                        phased_variants[(PS_flag_1,PS_flag_2)][hp_flag].append(data)
                        phased_variants_loci[(PS_flag_1,PS_flag_2)][hp_flag].append(int(data[1]))
                        #print(PS_flag_2 - PS_flag_1)

    count = 0
    count_2 = 0
    total_len = 0
    PS_dict = defaultdict(list)
    for key,val in phased_variants_loci.items():
        if len(val.keys()) == 2:
            count += 1
            min_loci_hp1 = min(val['hp1'])
            max_loci_hp1 = max(val['hp1'])
            min_loci_hp2 = min(val['hp2'])
            max_loci_hp2 = max(val['hp2'])
            min_loci = min([min_loci_hp1,min_loci_hp2])
            max_loci = max([max_loci_hp1,max_loci_hp2])
            PS_dict[key] = [min_loci,max_loci]
        else:
            count_2 += 1
            for key_2,val_2 in val.items():
                min_loci = min(val_2)
                max_loci = max(val_2)
                PS_dict[key] = [min_loci,max_loci]

    ### write to phased.vcf
    if chr_num == "X":
        chr_num = 23
    fw = open(out_dir + "phased_chr" + str(chr_num) + ".vcf","w")
    for key,val in phased_variants.items():
        print(key,key[1] - key[0])
        PS_info = PS_dict[key]
        PS_start = PS_info[0]
        PS_end = PS_info[1]
        for key_2, val_2 in val.items():
            if key_2 == "hp1":
                for data in val_2:
                    data[8] = "GT:PS"
                    data[9] = "1|0" + ":" + str(PS_start) + "-" + str(PS_end)
                    fw.write("\t".join(data) + "\n")
            elif key_2 == "hp2":
                for data in val_2:
                    data[8] = "GT:PS"
                    data[9] = "0|1" + ":" + str(PS_start) + "-" + str(PS_end)
                    fw.write("\t".join(data) + "\n")

    
    fw.close()
 

# 2. Phasing SNPs by Haplotyping Alg
def Phasing_SNPs_byAlg(in_dir,chr_num,SNP_dict,out_dir):
    fw = open(out_dir + "chr" + str(chr_num) + "_byAlg.vcf","w")
    # phasing blocks from haplotyping algorithms
    all_files = glob.glob(in_dir + "/" +  "chr" + str(chr_num) + "_*_final.p")
    all_pos = []
    total_num = 0
    for one_file in all_files:
        all_blocks = pickle.load(open(one_file,"rb"))
        total_num += len(all_blocks)
        for _num in range(len(all_blocks)):
            _block = all_blocks[_num]
            sorted_pos = sorted(_block.keys())
            _start_PS = np.min(sorted_pos)
            _end_PS = np.max(sorted_pos)
            print(_start_PS)
            for _pos, val in _block.items():
                all_pos.append(_pos)
                _key = (_pos + 1, val[1], val[2])
                if _key in SNP_dict:
                    data = SNP_dict[_key]
                    if val[0] == 2:
                        data[8] = data[8] + ":PS"
                        temp = data[9].split(":")
                        temp[0] = "1|0"
                        temp.append(str(_start_PS))
                        data[9] = ":".join(temp)
                        fw.write("\t".join(data) + "\n")
                    if val[0] == 1:
                        data[8] = data[8] + ":PS"
                        temp = data[9].split(":")
                        temp[0] = "0|1"
                        temp.append(str(_start_PS))
                        data[9] = ":".join(temp)
                        fw.write("\t".join(data) + "\n")
                    del SNP_dict[_key] 
                else:
                    print(_pos,val)

    fw.close()

    return len(all_pos)


def read_vcf_into_dict(vcf_file,chr_num,out_dir):
    if chr_num == 23:
        chr_num = "X"
    SNP_dict = defaultdict()
    with open(vcf_file) as f:
        for line in f:
            data = line.rsplit()
            if data[0][0] != "#":
                chr_num_curr = data[0]
                if chr_num_curr == "chr" + str(chr_num):
                    _pos = int(data[1])
                    ref = data[3]
                    alt = data[4]
                    GT = data[9].split(":")[0]
                    if (GT == "0/1" or GT == "1/0") and len(ref) == 1 and len(alt) == 1:
                    #if (GT == "0/1" or GT == "1/0"):
                        SNP_dict[(_pos,ref,alt)] = data
    if chr_num == "X":
        chr_num = 23
    pickle.dump(SNP_dict, open(out_dir + "chr" + str(chr_num) + "_SNPs_dict.p","wb"))
    return SNP_dict


# 3. Impute to phase assembled variants by SNPs(phased by Alg) 
def reverse_phase(phase_tag):
    if phase_tag == "0|1":
        return "1|0"
    elif phase_tag == "1|0":
        return "0|1"


def Compare_two_sets(phased_vcf,assembled_vcf,output_file,out_dir,chr_num,xin):
    temp_file = out_dir + "/_tmp_chr" + str(chr_num) + ".vcf"
    fw = open(temp_file,"w")
    dict_1 = defaultdict(list)
    dict_2 = defaultdict(list)
    with open(assembled_vcf,"r") as f:
        for line in f:
            data = line.rsplit()
            _pos = int(data[1]) + 1
            dict_1[_pos] = data
            
    with open(phased_vcf,"r") as f:
        for line in f:
            data = line.rsplit()
            _pos = int(data[1]) 
            dict_2[_pos] = data

    shared_loci = set(dict_1.keys()) & set(dict_2.keys())
    uniq_1 = list(set(dict_1.keys()) - shared_loci)
    uniq_2 = list(set(dict_2.keys()) - shared_loci)
    shared_loci_list = list(shared_loci)
    count = 0
    count_1 = 0
    count_wrong = 0
    count_miss_1 = 0
    count_miss_2 = 0
    for _pos in uniq_1:
        try:
            total_pos = shared_loci_list + [_pos]
            total_pos_sorted = sorted(total_pos)
            idx = total_pos_sorted.index(_pos)
            _pos_before = total_pos_sorted[idx - 1]
            _pos_after = total_pos_sorted[idx + 1]

            before_info_1 = dict_1[_pos_before][9].split(":")
            after_info_1 = dict_1[_pos_after][9].split(":")
            _info_1 = dict_1[_pos][9].split(":")
            
            before_info_2 = dict_2[_pos_before][9].split(":")
            after_info_2 = dict_2[_pos_after][9].split(":")

            if before_info_1[1] == after_info_1[1] and before_info_1[1] == _info_1[1] and before_info_2[-1] == after_info_2[-1]:
                before_phase_1 = before_info_1[0]
                after_phase_1 = after_info_1[0]
                _phase_1 = _info_1[0]
                before_phase_2 = before_info_2[0]
                after_phase_2 = after_info_2[0]
                if before_phase_1 == before_phase_2:
                    if after_phase_1 == after_phase_2:
                        count += 1
                        use_info = dict_1[_pos]
                        use_info[9] = _phase_1 + ":"  + before_info_2[-1]
                        use_info[1] =str(int(use_info[1]) + 1)   # + 1 for locus 
                        fw.write("\t".join(use_info) + "\n")
                    else:
                        count_wrong += 1
                if before_phase_1 != before_phase_2:
                    if after_phase_1 != after_phase_2:
                        count += 1
                        use_info = dict_1[_pos]
                        use_info[9] = reverse_phase(_phase_1) + ":"  + before_info_2[-1]
                        use_info[1] =str(int(use_info[1]) + 1)
                        fw.write("\t".join(use_info) + "\n")
                    else:
                        count_wrong += 1
            elif before_info_1[1] == _info_1[1]:
                count_1 += 1
                before_phase_1 = before_info_1[0]
                _phase_1 = _info_1[0]
                before_phase_2 = before_info_2[0]
                if before_phase_1 == before_phase_2:
                    use_info = dict_1[_pos]
                    use_info[9] = _phase_1 + ":"  + before_info_2[-1]
                    use_info[1] =str(int(use_info[1]) + 1)
                    fw.write("\t".join(use_info) + "\n")
                if before_phase_1 != before_phase_2:
                    use_info = dict_1[_pos]
                    use_info[9] = reverse_phase(_phase_1) + ":"  + before_info_2[-1]
                    use_info[1] =str(int(use_info[1]) + 1)
                    fw.write("\t".join(use_info) + "\n")
            elif after_info_1[1] == _info_1[1]:
                count_1 += 1
                after_phase_1 = after_info_1[0]
                _phase_1 = _info_1[0]
                after_phase_2 = after_info_2[0]
                if after_phase_1 == after_phase_2:
                    use_info = dict_1[_pos]
                    use_info[9] = _phase_1 + ":"  + after_info_2[-1]
                    use_info[1] =str(int(use_info[1]) + 1)
                    fw.write("\t".join(use_info) + "\n")
                if after_phase_1 != after_phase_2:
                    use_info = dict_1[_pos]
                    use_info[9] = reverse_phase(_phase_1) + ":"  + after_info_2[-1]
                    use_info[1] =str(int(use_info[1]) + 1)
                    fw.write("\t".join(use_info) + "\n")
            else:
                count_miss_1 += 1
        except:
            count_miss_2 += 1
            pass
    print("------------result for chr" + str(chr_num) + ": ")
    print(len(uniq_1))
    print(count,count_1,count_wrong)
    print(count_miss_1,count_miss_2)
    fw.close()
    correct_percent = float(count + count_1)/(count + count_1 + count_wrong)
    print(correct_percent)

    fw = open(output_file,"w")
    with open(temp_file,"r") as f:
        for line in f:
            data = line.rsplit()
            ref = data[3]
            alt = data[4]
            if ref == "-" and len(alt) == 1:
                data[7] = "SVTYPE=INS"
            elif alt == "-" and len(ref) == 1:
                data[7] = "SVTYPE=DEL"
            fw.write("\t".join(data) + "\n")

    fw.close()                


    return correct_percent


def get_homozygous_variants(vcf_file,output_file):
    fw = open(output_file,"w")
    with open(vcf_file,"r") as f:
        for line in f:
            data = line.rsplit()
            if data[0][0] != "#":
                GT = data[9].split(":")[0]
                if GT == "1|1" or GT == "1/1":
                    data[8] = "GT"
                    data[9] = "1|1"
                    fw.write("\t".join(data) + "\n")
    fw.close()




def main():
    if len(sys.argv) == 1:
        Popen("python3 " + "Aquila_phasing_all_variants.py -h",shell=True).wait()
    else:
        out_dir = args.out_dir + "/"
        if os.path.exists(out_dir):
            print("using existing output folder: " + out_dir)
        else:
            os.makedirs(out_dir)
        chr_start = args.chr_start
        chr_end = args.chr_end
        assembled_vcf = args.assembly_vcf
        vcf_file = args.vcf_file
        assemble_dir = args.assembly_dir
        phase_blocks_cut_highconf_dir = assemble_dir + "/phase_blocks_cut_highconf/"
        phased_file_dir = assemble_dir + "/results_phased_probmodel/"
        block_len_use = args.block_len_use
        # 1. phasing assembled variants
        for chr_num in range(chr_start,chr_end + 1):
            breakpoint_file = phase_blocks_cut_highconf_dir + "chr" + str(chr_num) + ".phased_final_cut_by_" + str(block_len_use) + "_HC_breakpoint_2.p" 
            HCbk_dict = pickle.load(open(breakpoint_file,"rb"))
            Phasing_assembled_SV(HCbk_dict,assembled_vcf,chr_num,out_dir)
        # 2. phasing SNPs by Alg
        total_num = 0
        total_snp = 0
        for chr_num in range(chr_start,chr_end + 1):
            SNP_dict = read_vcf_into_dict(vcf_file,chr_num,out_dir)
            ##SNP_dict = pickle.load(open(out_dir + "chr" + str(chr_num) + "_SNPs_dict.p","rb"))
            total_snp += len(SNP_dict)

            all_pos = Phasing_SNPs_byAlg(phased_file_dir, chr_num, SNP_dict,out_dir)
            total_num += all_pos
        print("total_num: " + str(total_num))
        print("total_snp: " + str(total_snp))
        # 3. Impute to phasing assembled variants
        pool = Pool(chr_end - chr_start + 1)
        for chr_num in range(chr_start,chr_end + 1):
            output_file = out_dir + "chr" + str(chr_num) + "_impute.vcf"
            phased_vcf = out_dir + "chr" +str(chr_num) + "_byAlg.vcf"
            assembled_vcf = out_dir + "phased_chr" + str(chr_num) + ".vcf"
            # pool.apply_async(Compare_two_sets,(phased_vcf,assembled_vcf,output_file,out_dir,chr_num,"xin"))
            Compare_two_sets(phased_vcf,assembled_vcf,output_file,out_dir, chr_num,"xin")
        pool.close()
        while len(active_children()) > 1:
            time.sleep(0.5)
        pool.join()
        print("done~")

        for chr_num in range(chr_start,chr_end + 1):
            one_impute = out_dir + "chr" + str(chr_num) + "_impute.vcf"
            one_byAlg = out_dir + "chr" + str(chr_num) + "_byAlg.vcf"
            cat_cmd_1 = "cat " + one_impute  + " >> " + out_dir + "all_impute.vcf"
            cat_cmd_2 = "cat " + one_byAlg  + " >> " + out_dir + "all_byAlg.vcf"
            Popen(cat_cmd_1,shell=True).wait()
            Popen(cat_cmd_2,shell=True).wait()

        homo_file = out_dir + "all_homo_assembled.vcf"
        get_homozygous_variants(assembled_vcf,homo_file)

        final_vcf = out_dir + "Aquila_all_phased.vcf"
        final_sorted_vcf = out_dir + "Aquila_all_phased_sorted.vcf"
        final_cat_cmd = "cat " + out_dir + "all_impute.vcf " + out_dir + "all_byAlg.vcf " + homo_file + " > " + final_vcf
        final_sort_cmd = "cat " + final_vcf + " | awk '$1 ~ /^#/ {print $0;next} {print $0 | \"LC_ALL=C sort -k1,1 -k2,2n\"}' > " + final_sorted_vcf
        Popen(final_cat_cmd,shell=True).wait()
        Popen(final_sort_cmd,shell=True).wait()


if __name__ == "__main__":
    main()
