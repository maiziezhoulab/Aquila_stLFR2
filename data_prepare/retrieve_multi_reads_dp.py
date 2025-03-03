from collections import defaultdict 
idxl = [11, 15, 18, 19, 21, 22, 23, 24, 25, 26, 27]

dc = defaultdict(dict)

with open("/data/maiziezhou_lab/ZhouLab_Projects/VandyCG_stLFR/Experiments/2023_October/Aquila_stLFR_chr6/read/unused_name_uniq_sorted.txt",'r') as f:
    for line in f:
        rnv = [line[:-1][i]  for i in idxl]


        last_dc = defaultdict(list)
        last_dc[rnv[-2]].append(rnv[-1])
        print(last_dc)
        
        print(len(rnv))
        for i in range(len(rnv)-2,1,-1):
            new_dc = defaultdict(dict)
            new_dc[rnv[i-1]] = last_dc
            last_dc = new_dc.copy()
            print(i,last_dc)
        dc[rnv[0]] = last_dc
        break 

print(rnv)
print(dc)
l = ['1', '1', '0', '1', '0','0', '0', '0', '6', '4', '5']
print(l)
dc1 = dc.copy()
for x in l[:-1]:
    dc1 = dc1[x].copy()
    print(x,dc1)
# cnt = 0
# with open("/data/maiziezhou_lab/ZhouLab_Projects/VandyCG_stLFR/Experiments/2023_October/Aquila_stLFR_chr6/read/unused_name_uniq_sorted.txt",'r') as f:
#     for line in f:
#         cnt+=1
#         if cnt<=1:
#             continue

#         rnv = [line[:-1][i]  for i in idxl]
#         for r in rnv:





