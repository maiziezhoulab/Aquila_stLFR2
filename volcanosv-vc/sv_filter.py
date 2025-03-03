import pysam 
import argparse
import numpy as np
from argparse import ArgumentParser
# from subprocess import Popen
import subprocess
import os
parser = ArgumentParser(description="",
	usage='use "python3 %(prog)s --help" for more information',
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--vcffile','-v')
parser.add_argument('--bamfile','-b')
parser.add_argument('--bamcov','-cov', )
# parser.add_argument('--output_path','-o')
# parser.add_argument('--dtype','-d', choices= ['SE','PE'])
parser.add_argument('--n_thread','-t', type = int, default = 50 )
parser.add_argument('--ratio','-r', type = float, default = 0.84)
# parser.add_argument('--delete_temp_file','-d', action='store_true')
args = parser.parse_args()
vcffile = args.vcffile
bamfile = args.bamfile
bamcov = eval(args.bamcov)
# output_path = args.output_path
# dtype = args.dtype
n_thread = args.n_thread
r = args.ratio
# chr_num = args.chr_num


thresh = bamcov * r
print("thresh ", thresh)
output_path = os.path.abspath(vcffile).replace(".vcf",'_filtered.vcf')
def load_vcf(vcffile):
    header = []
    dels = []
    with open(vcffile,'r') as f:
        for line in f:
            if line[0]=='#':
                header.append(line)
            elif 'SVTYPE=DEL' in line:
                data = line.split()
                svlen = len(data[4]) - len(data[3])
                chrom = data[0]
                pos = int(data[1])
                end = pos + abs(svlen)
                svid = data[2]
                var = (chrom, pos, end, svid, line)
                dels.append(var)
            else:
                header.append(line)
    return dels, header

def check_coverage_se(bamfile,var, flanking):
    chrom,start_pos,end_pos,svid, _ = var
    svlen = end_pos-start_pos
    samfile = pysam.AlignmentFile(bamfile)
    # dc = defaultdict(list)
    cov_list = []
    dc = {}
    for pileupcolumn in samfile.pileup(chrom, start_pos - flanking, end_pos+ flanking, min_base_quality = 0, truncate = True): 
        cov = 0
        for pileupread in pileupcolumn.pileups:
            if pileupread.is_del or pileupread.is_refskip:
                # Deletion or reference skip

                continue
            else:
                cov +=1
        cov_list.append(cov)
        dc[pileupcolumn.pos] = cov 


    for pos in range(start_pos - flanking, end_pos + flanking):
        if pos  not in dc:
            dc[pos] = 0 
    return dc,svid,start_pos,end_pos

def check_coverage(bamfile,var, flanking):
    chrom,start_pos,end_pos,svid, varline = var
    svlen = end_pos-start_pos

    result = subprocess.run(['samtools', 'depth','-a', bamfile,'-r', f"{chrom}:{start_pos-flanking}-{end_pos+flanking}"], stdout=subprocess.PIPE)
    out = result.stdout.decode('utf-8')
    # print(out)
    # exit()

    dc = {}
    for line in out.split('\n')[:-1]:
        _, pos, cov = line.split()
        cov = int(cov)
        pos = int(pos)
        dc[pos] = cov
    

    return dc,svid,start_pos,end_pos,varline


dels, header = load_vcf(vcffile)
print(len(dels))
flanking = 30
from joblib import Parallel, delayed
from tqdm import tqdm
# n_thread = 10
# assert dtype in ['PE','SE']
# if dtype =='PE':
#     sequences = Parallel(n_jobs=n_thread)(delayed(check_coverage_pe )(bamfile,var, flanking) 
#                                         for var in tqdm(dels))
# else:
#     sequences = Parallel(n_jobs=n_thread)(delayed(check_coverage_se )(bamfile,var, flanking) 
#                                         for var in tqdm(dels))
   
sequences = Parallel(n_jobs=n_thread)(delayed(check_coverage )(bamfile,var, flanking) 
                                        for var in tqdm(dels))

ft_cnt = 0
with open(output_path,'w') as f:
    f.writelines(header)
    for dc,svid,start_pos,end_pos,varline in sequences:
        s = f'{svid}\t{start_pos}\t{end_pos}'
        cov_list = []
        for pos in range(start_pos- flanking, end_pos + flanking):
            if pos in dc:
                cov = dc[pos]
            else:
                cov = 0
            x = f'\t{pos}:{cov}'
            s = s+x 
            cov_list.append(cov)

        cov_list = np.array(cov_list)
        delta = sorted(abs(cov_list[:-1] - cov_list[1:]))
        peaks = np.array(delta[-2:]).mean()
        mean_cov = cov_list[flanking:-flanking].mean()
        # print(mean_cov)
        # if mean_cov<= 40:
        if mean_cov <= thresh:
            f.write(varline)
        
        if peaks<=2:
            ft_cnt+=1

        s = s +'\n'

        # s = svid+','+str(svlen)+';'+','.join([str(x) for x in seq])+'\n'
        # print(seq)
        # f.write(s)


# print(ft_cnt)

# cmd = f'''/data/maiziezhou_lab/CanLuo/long_reads_project/bin/truvari_eval_v1.sh {chr_num} \
#     {os.path.dirname(output_path)}  \
#         {os.path.dirname(output_path)}/eval_filtered/ \
#               {os.path.basename(output_path).split('.')[0]} \
#               500 0.5 0.5 30 0.01
#         '''
# print(cmd)
# os.system(cmd)


# cmd = f'''/data/maiziezhou_lab/CanLuo/long_reads_project/bin/truvari_eval_v1.sh {chr_num} \
#     {os.path.dirname(output_path)}  \
#         {os.path.dirname(output_path)}/eval_filtered/ \
#               {os.path.basename(output_path).split('.')[0]} \
#               500 0 0 30 0
#         '''
# os.system(cmd)