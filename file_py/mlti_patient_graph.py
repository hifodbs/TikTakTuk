# %%
import pyreadr

# Load the rds file
result = pyreadr.read_r('06_Cb_BTM_table.rds')

# pyreadr returns a dictionary; the data is usually under the 'None' key 
# because rds files contain a single object without a specific name.
df = result[None]

# select ttype 

ttype = "OV"

df = df[df["ttype"]==ttype]

df = df.dropna(subset=["clock_rank","gene"])

trees = {}

patients = df["sample_id"].unique()

interesting_p = 0

for p in patients:
    df_p = df[df["sample_id"] == p]
    ranks = int(df_p["clock_rank"].unique().max())
    genes = {}
    for r in range(1,ranks+1):
        df_p_r =  df_p[df_p["clock_rank"]==r]
        genes_list = df_p_r["gene"].tolist()
        if len(genes_list)>0: #some ranks are skipped ????
            genes[r] = df_p_r["gene"].tolist()
    trees[p] = genes
    interesting_p += 1
    

# %%
# tring some mixup
from collections import defaultdict
from utils import plot_tree as pt
import math


# map gene end in which clock was found with multeplicity

r_gen_all = defaultdict(lambda: defaultdict(int))

for p, ranks in trees.items():
    for r, genes_list in ranks.items():
        for g in genes_list:
            r_gen_all[g][r] += 1
            
print("All gens found : ",len(r_gen_all))


# take only the rank where it was found most of the time and remap rank -> (gene, multiplicity)
r_gen_refined = defaultdict(list)
for k,v in r_gen_all.items():
    r_gen_refined[max(v, key=v.get)].append((k,v[max(v, key=v.get)]))
    
    
# drop genes that appear only less than treshold

temp = {}
for r,genes_list in r_gen_refined.items():
    genes_list.sort(key=lambda tup: tup[1], reverse=True)
    if len(genes_list)>20:
        genes_list = genes_list[0:int(math.sqrt(len(genes_list)))]
    else:
        genes_list = genes_list[0:5]
    temp[r] = genes_list
r_gen_refined = temp
        
print(r_gen_refined[1])
        
    
multi_tree = defaultdict(list)
for r,genes_list in r_gen_refined.items():
    for gene, mult in genes_list:
        multi_tree[r].append(gene)


    
pt.print_tree(multi_tree,"multi_patient.gv")


