
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


def extract_vcf_info(vcf_file,chrom, out_dir, mbq_threshold):
    f = open(vcf_file,"r")
    curr = 0
    variant_dict =  defaultdict(int)
    variant_dp_dict =  defaultdict(int)
    curr = 0
    count_total = 0
    total_depth = []
    use_chrs = set([f'chr{i}' for i in range(1,23)])
    for line in f:
        data = line.rsplit()
        curr += 1
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





def split_chrom(bam, chrom, n):
    b = pysam.AlignmentFile(bam, "rb")
    L = b.get_reference_length(chrom)
    s = math.ceil(L / n)
    return [(i*s, min((i+1)*s, L)) for i in range(n)]


def segment_sites_by_ranges(sites, ranges):
    result = []
    for start, end in ranges:
        block = [(p, b) for p, b in sites if start <= p < end]
        result.append(block)
    return result


def split_pos_ref_into_blocks(pos_ref_list, n_blocks):
    """
    Split a list of (position, ref_allele) tuples into n_blocks such that
    each block covers roughly equal reference span (max(pos)-min(pos)).

    Parameters
    ----------
    pos_ref_list : list[tuple[int, str]]
        List of (position, ref_allele)
    n_blocks : int
        Number of desired blocks

    Returns
    -------
    list[list[tuple[int, str]]]
        List of n_blocks (or fewer if too few positions)
    """
    if not pos_ref_list:
        return []
    if n_blocks <= 1:
        return [sorted(pos_ref_list, key=lambda x: x[0])]

    # Sort by position
    pos_ref_list = sorted(pos_ref_list, key=lambda x: x[0])
    positions = [p for p, _ in pos_ref_list]

    total_span = positions[-1] - positions[0]
    target_span = total_span / n_blocks

    blocks = []
    current_block = [pos_ref_list[0]]
    block_start = positions[0]

    for (pos, ref) in pos_ref_list[1:]:
        if pos - block_start <= target_span:
            current_block.append((pos, ref))
        else:
            blocks.append(current_block)
            current_block = [(pos, ref)]
            block_start = pos

    if current_block:
        blocks.append(current_block)

    # Merge small blocks if we created too many
    while len(blocks) > n_blocks:
        spans = [b[-1][0] - b[0][0] for b in blocks]
        i = spans.index(min(spans))
        if i < len(blocks) - 1:
            blocks[i].extend(blocks[i+1])
            del blocks[i+1]
        else:
            blocks[i-1].extend(blocks[i])
            del blocks[i]

    return blocks

def per_read_alleles_refpos(bam_path, chrom, block, sites, out_dir, task_id,  min_mapq=20, one_based=False, maxd=50_000):
    """
    Extracts per-read alleles at specific reference positions using get_reference_positions().

    Parameters:
        bam_path : str
            Path to BAM file.
        sites : list of (chrom, pos)
            Genomic sites (0-based by default; use one_based=True if your list is 1-based).
        min_mapq : int
            Minimum mapping quality.
        one_based : bool
            Whether positions are 1-based.

    Returns:
        pandas.DataFrame with columns:
        [chrom, pos, read, allele, is_del, read_start, strand, mapq, bx]
    """
    bam = pysam.AlignmentFile(bam_path, "rb")
    rows = []


    pos_ref_mapper = dict(sites)
    pos_set = set([ var[0] for var in sites])

    mole_read_dict = defaultdict(list)


    for read in bam.fetch(chrom, block[0], block[1] + 1):
        if read.is_unmapped or read.mapping_quality < min_mapq:
            continue

        ref_positions = read.get_reference_positions(full_length=True)
        seq = read.query_sequence
        # assert len(seq) == len(ref_positions)
        if len(ref_positions) != len(read.query_sequence):
            print(read.query_name, read.pos, read.cigarstring, len(read.query_sequence), len(ref_positions))
        bx = read.get_tag("BX") if read.has_tag("BX") else None

        which_read = 1 if read.is_read1 else 2


        cur_read_var = []
        for i, refpos in enumerate(ref_positions):
            if refpos in pos_set:
                # None means insertion or soft clip, skip
                if refpos is None or i >= len(seq):
                    continue
                base = seq[i]
                ref_allele = pos_ref_mapper[refpos]
                if base == ref_allele:
                    allele = 0
                else:
                    allele = 1

                cur_snp_pos =  refpos + 1 if one_based else refpos
                cur_read_var.append((cur_snp_pos,allele))

        if bx not in mole_read_dict:

            mole_read_dict[bx] = [[(read.qname, which_read, read.pos, cur_read_var)]]
        else:
            last_pos = mole_read_dict[bx][-1][-1][2]
            dist = read.pos - last_pos
            if dist > maxd:

                mole_read_dict[bx].append([(read.qname, which_read, read.pos, cur_read_var)])
            else:

                mole_read_dict[bx][-1].append((read.qname, which_read, read.pos, cur_read_var))
    bam.close()
    outfile = f"{out_dir}/mole_dict_{task_id}_{chrom}_{block[0]}_{block[1]}.p"
    with open(outfile,'wb') as f:
        pickle.dump(mole_read_dict, f)
    return 

def organize_snps(read_list):
    dc = defaultdict(list)
    for read in read_list:
        cur_var = read[3]
        for pos, allele in cur_var:
            dc[pos].append(allele)
    return dc 

def check_hp_one_snp(hp,):
    if len(set(hp)) == 1:
        hp_use = True
        hp_flag = hp[0]
    else:
        count0 = hp.count(0)
        count1 = hp.count(1)
        if count0 >= 4*count1:
            hp_use = True
            hp_flag = "0"
        elif count1 >= 4*count0:
            hp_use = True
            hp_flag = "1"
        else:
            hp_use = False
            hp_flag = "none"
    return (hp_use,hp_flag)

def check_hp(mole_variant):
    s = []
    for locus_start, hp in mole_variant.items():
        hp_use, hp_flag = check_hp_one_snp(hp)
        if hp_use:
            s.append(f"{locus_start}:{hp_flag}")
    if len(s):
        s = '\t'.join(s)+'\n'
    else:
        s = '\n'
    return s

def parse_mole(mole, read_type, chrom, bx,mole_id):
    dc = defaultdict(list)
    dc_qname_pos = defaultdict(list)
    for read in mole:
        dc[read[0]].append(read[1])
        dc_qname_pos[read[0]].append(read[2])
    if read_type == 'PE':
        valid_pair = 0 
        for k,v in dc.items():
            if len(set(v))==2:
                valid_pair+=1 
    else:
        valid_pair = len(dc)

    if valid_pair>=2 :
        # dc_read = defaultdict(list)
        mole_variant = organize_snps(mole)
        hp_str = check_hp(mole_variant)

        # for read in mole_variant:
        #     dc_read[read[0]].append(read[1])
        poss = [read[2] for read in mole]
        line = f"{chrom[3:]}\t{min(poss)+1}\t{max(poss)+1}\t{max(poss)-min(poss)+100}\t{len(poss)}\t{bx}\t{mole_id}\t{hp_str}"
    else:
        line = ""
        # dc_read = defaultdict(list)/
    return line, dc, dc_qname_pos
        
def process_mole_dict(mole_dict, read_type, chrom, start_id, fw):
    mole_id = start_id
    mole_qname_dict = {}
    qname_pos = defaultdict(list)
    for bx, mole_list in mole_dict.items():
        for mole in mole_list:
            line, dc_read, dc_qname_pos = parse_mole(mole, read_type, chrom, bx, mole_id)
            if line!='':
                mole_qname_dict[mole_id] = dc_read

                for qname, poss in dc_qname_pos.items():
                    qname_pos[qname].extend(poss)
                # print(line)
                mole_id+=1 
                fw.write(line)

    return mole_id, mole_qname_dict, qname_pos

def process_one_pkl(pkl_list, i,  chrom, maxd):

    if i==0:
        prev_pkl = None 
        cur_pkl = pkl_list[i]
        next_pkl = pkl_list[i+1]
    elif i == (len(pkl_list)-1):
        prev_pkl = pkl_list[i-1] 
        cur_pkl = pkl_list[i]
        next_pkl = None 
    else:
        prev_pkl = pkl_list[i-1] 
        cur_pkl = pkl_list[i]
        next_pkl = pkl_list[i+1]

        
    if prev_pkl is not None:
 
        d1 = pickle.load(open(prev_pkl, 'rb'))
        d2 = pickle.load(open(cur_pkl, 'rb'))
        if (len(d1) == 0) | (len(d2)==0):
            cur_dc = d2 
        else:
            k1, k2 = set(d1), set(d2)
            inter = k1 & k2
            print(f"A:{len(k1)} B:{len(k2)} overlap:{len(inter)} "
                f"A_in_B:{len(inter)/len(k1)*100:.2f}% B_in_A:{len(inter)/len(k2)*100:.2f}%")
            cnt = 0
            for bx in inter:
                prev_end = d1[bx][-1][-1][2]
                cur_start = d2[bx][0][0][2]
                if (cur_start - prev_end) < maxd:
                    prev_len = len(d1[bx])
                    d2[bx][0] = d1[bx][-1] + d2[bx][0]
                    d1[bx].pop()
                    cur_len = len(d1[bx])
                    assert cur_len ==  (prev_len -1)
                    cnt+=1

            cur_dc = d2
            print(cnt)
    else:
        cur_dc = pickle.load(open(cur_pkl, 'rb'))

    if next_pkl is not None:
        d1 = cur_dc 
        d2 = pickle.load(open(next_pkl, 'rb'))
        if (len(d1) == 0) | (len(d2)==0):
            cur_dc = d1 
        else:
            k1, k2 = set(d1), set(d2)
            inter = k1 & k2
            print(f"A:{len(k1)} B:{len(k2)} overlap:{len(inter)} "
                f"A_in_B:{len(inter)/len(k1)*100:.2f}% B_in_A:{len(inter)/len(k2)*100:.2f}%")
            cnt = 0
            for bx in inter:
                prev_end = d1[bx][-1][-1][2]
                cur_start = d2[bx][0][0][2]
                if (cur_start - prev_end) < maxd:
                    prev_len = len(d1[bx])
                    d2[bx][0] = d1[bx][-1] + d2[bx][0]
                    d1[bx].pop()
                    cur_len = len(d1[bx])
                    assert cur_len ==  (prev_len -1)
                    cnt+=1
            print(cnt)
            cur_dc = d1

    outfile = cur_pkl+".h5"
    with open(outfile, 'w') as fw:
        process_mole_dict(cur_dc, read_type, chrom, 0, fw)
    
def safe_update(qname_dict, qname_dict_cur):

    olp_qnames = set(list(qname_dict_cur.keys())) & set(list(qname_dict.keys()))

    for qname in olp_qnames:
        qname_dict[qname].extend(qname_dict_cur[qname])
        del qname_dict_cur[qname]

    qname_dict.update(qname_dict_cur)
    return qname_dict


def write_h5(pkl_list , out_dir, sample, chrom,  read_type, maxd= 50_000):

    start_id = 0

    h5_file = f"{out_dir}/{sample}_{chrom}.h5"
    mole_qname_file = f"{out_dir}/{sample}_{chrom}_qname.p"
    qname_pos_file = f"{out_dir}/{sample}_{chrom}_qname_pos.p"

    fw = open(h5_file,'w')

    cur_dict = pickle.load(open(pkl_list[0], 'rb'))

    for i in tqdm(range(len(pkl_list)-1), desc = "write molecule to h5"):
        pkl1 = pkl_list[i]
        pkl2 = pkl_list[i+1]

        d1 = cur_dict
        # d1 = pickle.load(open(pkl1, 'rb'))
        d2 = pickle.load(open(pkl2, 'rb'))
        if (len(d1)==0) or (len(d2)==0):
            pass
        else:

            k1, k2 = set(d1), set(d2)
            inter = k1 & k2
            print(f"A:{len(k1)} B:{len(k2)} overlap:{len(inter)} "
                f"A_in_B:{len(inter)/len(k1)*100:.2f}% B_in_A:{len(inter)/len(k2)*100:.2f}%")
            cnt = 0
            for bx in inter:
                prev_end = d1[bx][-1][-1][2]
                cur_start = d2[bx][0][0][2]
                if (cur_start - prev_end) < maxd:
                    prev_len = len(d1[bx])
                    d2[bx][0] = d1[bx][-1] + d2[bx][0]
                    d1[bx].pop()
                    cur_len = len(d1[bx])
                    assert cur_len ==  (prev_len -1)
                    cnt+=1
            print(cnt)
        
        end_id, mole_qname_dict_batch, qname_pos_batch = process_mole_dict(d1, read_type, chrom, start_id, fw)
        if i == 0:
            mole_qname_dict = mole_qname_dict_batch
            qname_pos = qname_pos_batch 
        else:
            mole_qname_dict.update(mole_qname_dict_batch)
            # qname_pos.update(qname_pos_batch)
            qname_pos = safe_update(qname_pos, qname_pos_batch)
        start_id = end_id+1

        cur_dict = d2 

    # with open(pkl_list[-1]+".new",'wb') as f:
    #     pickle.dump(d2,f)

    end_id, mole_qname_dict_batch, qname_pos_batch = process_mole_dict(d2, read_type, chrom, start_id, fw)
    mole_qname_dict.update(mole_qname_dict_batch)
    # qname_pos.update(qname_pos_batch)
    qname_pos = safe_update(qname_pos, qname_pos_batch)
    fw.close()

    with open(mole_qname_file,'wb') as f:
        pickle.dump(mole_qname_dict, f)

    with open(qname_pos_file,'wb') as f:
        pickle.dump(qname_pos, f)

def merge(outdir, outfile):
    cmd = '''cat %s/mole_dict*.p.h5 | sort -k2,2n | awk 'BEGIN{FS=OFS="\t"} {$7=NR; print}' > %s'''%(outdir, outfile)
    print(cmd)
    Popen(cmd, shell= True).wait()


def sort(infile):
    outfile = infile.replace('.h5','_sorted.h5')
    cmd = '''cat %s | sort -k2,2n  > %s'''%(infile, outfile)
    print(cmd)
    Popen(cmd, shell= True).wait()


def gen_H5(bamfile, vcf_file, sample, chr_num,  read_type, outdir, n_thread, mbq_threshold = 13, min_mapq = 20, max_d =50_000):
    
    n_blocks = n_thread 
    chrom = f'chr{chr_num}'

    start = time.time()

    if not os.path.exists(outdir):
        os.system("mkdir -p " + outdir)
    extract_vcf_info(vcf_file,chrom, outdir, mbq_threshold)

    dc = pickle.load(open(f"{outdir}/variant_dict_heterozygous.p",'rb'))

    all_vars = []
    for key,val in dc.items():
        if key[0]==chr_num:
            all_vars.append((key[1], val[0]))


    print(len(dc), len(all_vars))

    blocks = split_chrom(bamfile, chrom, n_blocks)

    var_blocks = segment_sites_by_ranges(all_vars, blocks)

    pkl_list = [f"{outdir}/mole_dict_{i}_{chrom}_{block[0]}_{block[1]}.p"for i,block in enumerate(blocks)]


    results = Parallel(n_jobs=n_thread)(delayed(per_read_alleles_refpos)\
                                        (bamfile, chrom, blocks[task_id], var_blocks[task_id], outdir, task_id, min_mapq, False, max_d) 
                                        for task_id in tqdm(range(len(blocks))))

    write_h5(pkl_list, outdir, sample, chrom, read_type )
    sort(f"{outdir}/{sample}_{chrom}.h5")

    end = time.time()

    print(f"Elapsed time: {end - start:.2f} seconds")



# n_blocks = 40

# outdir = 'temp'
# read_type = 'PE'
# sample = 'HG002'
# chr_num = 1
# vcf_file = "/lfs/archer.accre.vu/maiziezhou_lab/maiziezhou_lab/Datasets/CG_datasets/V350158928/BWA_align/gatk/gatk.vcf"
# bamfile = "/lfs/archer.accre.vu/maiziezhou_lab/maiziezhou_lab/Datasets/CG_datasets/V350158928/BWA_align/V350158928_BWA_hg19.bam"


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="",
        usage='use "python3 %(prog)s --help" for more information',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser = argparse.ArgumentParser(description="Extract reads for all small chunks:")
    parser.add_argument('--bam_file','-bam',help="bam file")
    parser.add_argument('--vcf_file','-vcf',help="vcf file")
    parser.add_argument('--sample_name','-s',help="sample name")
    parser.add_argument('--out_dir','-o',help="output folder")
    parser.add_argument('--read_type','-r',help = 'SE/PE', choices= ['SE','PE'])
    parser.add_argument('--mbq_threshold','-bq',type=int,help="mbq threshold", default=13)
    parser.add_argument('--mmq_threshold','-mq',type=int,help="mmq threshold", default=20)
    parser.add_argument('--chr_num','-chr',type=int,help="chromosome number", default=1)
    parser.add_argument('--boundary','-b',type=int,help="cut boundary", default=50_000)
    parser.add_argument('--num_threads','-n',type=int,help="number of threads", default=8)
    

    args = parser.parse_args()

    bam_file = args.bam_file
    vcf_file = args.vcf_file
    mbq_threshold = args.mbq_threshold
    mmq_threshold = args.mmq_threshold
    sample_name = args.sample_name
    chr_num = args.chr_num
    num_threads = args.num_threads
    read_type = args.read_type
    out_dir = args.out_dir
    boundary = args.boundary
    assert read_type in ['SE','PE']

    gen_H5(bam_file, vcf_file, sample_name, chr_num,  read_type, out_dir, num_threads, mbq_threshold , mmq_threshold , boundary )

