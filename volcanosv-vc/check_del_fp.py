import argparse
from argparse import ArgumentParser
parser = ArgumentParser(description="",
	usage='use "python3 %(prog)s --help" for more information',
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--vcffile','-v')
parser.add_argument('--test_range','-r', help = "a-b")
parser.add_argument('--del_eval_dir','-e')
parser.add_argument('--outfile','-o', help = 'csv')
# parser.add_argument('--n_thread','-t', type = int, default = 22 )
# parser.add_argument('--delete_temp_file','-d', action='store_true')
args = parser.parse_args()
vcffile = args.vcffile
test_range = args.test_range.split('-')
del_eval_dir = args.del_eval_dir
outfile = args.outfile
# n_thread = args.n_thread

start = int(test_range[0])
end = int(test_range[1])
import pandas as pd
import numpy as np
from tqdm import tqdm
def load_vcf(vcffile):
    dc ={}
    with open(vcffile,'r') as f:
        for line in f:
            if (line[0]!='#') & ('DEL' in line):
                data = line.split()
                svid = data[2]
                cov = eval(data[-1])
                dc[svid] = cov
    return dc 

def load_sum(sumfile):
    with open(sumfile,'r') as f:
        dc = eval(f.read())
    return dc 


def load_vcf_svid(vcffile):
    svid_list = []
    with open(vcffile,'r') as f:
        for line in f:
            if line[0]!='#':
                svid = line.split()[2]
                svid_list.append(svid)

    return svid_list 



tp_file = del_eval_dir+'/tp-call.vcf'
fp_file = del_eval_dir+'/fp.vcf'
sum_file = del_eval_dir+'/summary.txt'

tp_sv = load_vcf_svid(tp_file)
fp_sv = load_vcf_svid(fp_file)
dc_sum = load_sum(sum_file)
dc_sv = load_vcf(vcffile)
print(dc_sum)
data = []
for svid,cov in dc_sv.items():
    if svid in set(tp_sv):
        data.append([svid,cov,'tp'])
    elif svid in set(fp_sv):
        data.append([svid,cov,'fp'])

df = pd.DataFrame(data, columns = ['svid','cov','type']).sort_values('cov').reset_index(drop = True)
cov_list = df['cov'].values
type_list = df['type'].values 

total_bench = dc_sum['TP-base'] + dc_sum['FN']
total_call = dc_sum['TP-call'] + dc_sum['FP']
assert total_call == df.shape[0]

data = []
for i in tqdm(range(start,end)):
    keep = np.where(cov_list<=i)[0]
    keep_list = type_list[keep]
    new_tp = (keep_list == 'tp').sum()
    new_fp = (keep_list == 'fp').sum()
    new_total_call = new_tp+new_fp
    new_recall = np.round(new_tp / total_bench,2)
    new_precision = np.round(new_tp/new_total_call,2)
    new_f1 = np.round(2 * new_precision * new_recall / (new_precision + new_recall),2)
    data.append([i,new_tp,new_fp,new_total_call, total_bench, new_recall, new_precision, new_f1])

df = pd.DataFrame(data, columns = ['thresh','new_tp','new_fp',
'new_total_call', 'total_bench', 'new_recall', 'new_precision', 'new_f1'])


df = df.sort_values(['new_f1','new_recall'],ascending = False).reset_index(drop = True)

print(df.head())

df.to_csv(outfile, index = False)




















