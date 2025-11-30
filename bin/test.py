import pandas as pd
from collections import defaultdict
import pickle
import numpy as np
import os
from subprocess import Popen
import time
from joblib import Parallel, delayed
from tqdm import tqdm
import pickle
import pysam, math
import psutil
import gc
def extract_vcf_info(vcf_file,chrom, out_dir, mbq_threshold):
    f = open(vcf_file,"r")
    curr = 0
    variant_dict =  defaultdict(int)
    variant_dp_dict =  defaultdict(int)
    curr = 0
    count_total = 0
    total_depth = []
    use_chrs = set([f'chr{i}' for i in range(1,23)])
    try:
        for line in f:
            data = line.rsplit()
            curr += 1
            # print(curr)
            if data[0][:2] != "##" and data[0][:2] != "#C" and data[0] in use_chrs:
                chr_num = int(data[0][3:])
                pos = int(data[1]) - 1
                ref = data[3]
                alt = data[4]
                qual = float(data[5])
                GT = data[9].split(":")[0]
                if (GT == "0/1"  or GT == "1/0") and len(ref) == 1 and len(alt) == 1 and qual >= mbq_threshold:
                    if data[0] == chrom:
                        variant_dict[(chr_num,pos)] = [ref,alt,GT]
                    break_flag = 0
                    _format = data[8].split(":")
                    DP_idx = _format.index("DP")
                    _format_info = data[9].split(":")
                    dp_depth = float(_format_info[DP_idx])

                    if 'AO' in _format:
                        # freebayes format
                        AO_idx = _format.index("AO")
                        RO_idx = _format.index("RO")
                        ao_depth = float(_format_info[AO_idx])
                        ro_depth = float(_format_info[RO_idx])
                    else:
                        # gatk format
                        AD_idx = _format.index('AD')
                        AD = _format_info[AD_idx]
                        ro_depth, ao_depth = AD.split(',')
                        ro_depth = float(ro_depth)
                        ao_depth = float(ao_depth)
                    
                    if ao_depth >= ro_depth:
                        ratio = ro_depth/ao_depth
                    else:
                        ratio = ao_depth/ro_depth

                    for ii in range(0,11):
                        if dp_depth <= 10*(ii + 1) and dp_depth > 10*ii: 
                            for jj in range(0,20):
                                if ratio <= 0.05*(jj+1) and ratio > 0.05*jj:
                                    variant_dp_dict[(ii,jj)] += 1
                                    count_total += 1
                                    break_flag = 1
                                    total_depth.append(dp_depth)
                                    break
                    if break_flag == 0:
                        if dp_depth > 110:
                            for jj in range(0,20):
                                if ratio <= 0.05*(jj+1) and ratio > 0.05*jj:
                                    variant_dp_dict[(11,jj)] += 1
                                    count_total += 1
                                    total_depth.append(dp_depth)
    except:
        print(curr)
    f.close()
    pickle.dump(variant_dict, open(out_dir + "/variant_dict_heterozygous.p", "wb"))
    avg_depth = np.mean(total_depth)
    median_depth = np.median(total_depth)
    print(count_total)
    print(np.mean(total_depth))
    print(np.median(total_depth))
    print(np.max(total_depth))
    print(np.min(total_depth))

    var_depth_file = out_dir + "/median_depth_for_var.txt"
    f = open(var_depth_file,"w")
    f.writelines(str(median_depth) + "\n")
    f.close()
    return 

# vcf_file = "/lfs/archer.accre.vu/maiziezhou_lab/maiziezhou_lab/Datasets/CG_datasets/V350158968/BWA_align/gatk/gatk.vcf"
vcf_file = "/lfs/archer.accre.vu/maiziezhou_lab/maiziezhou_lab/Datasets/CG_datasets/V350158968/BWA_align/gatk/new"
# vcf_file = "/lfs/archer.accre.vu/maiziezhou_lab/maiziezhou_lab/Datasets/CG_datasets/V350159389/BWA_align/gatk/gatk.vcf"
chrom, out_dir, mbq_threshold = 'chr1','out',20
extract_vcf_info(vcf_file,chrom, out_dir, mbq_threshold)