# %%
import pyreadr

# Load the rds file
result = pyreadr.read_r('06_Cb_BTM_table.rds')

# pyreadr returns a dictionary; the data is usually under the 'None' key 
# because rds files contain a single object without a specific name.
df = result[None]

df = df[df['mutatation_status']!="WT"]

# select ttype 

ttype = "ESAD"

df = df[df["ttype"]==ttype]

df = df.dropna(subset=["clock_rank","gene"])

df

# %%
trees = {}

patients = df["sample_id"].unique()

interesting_p = 0
max_gene = 0


for p in patients:
    df_p = df[df["sample_id"] == p]
    genes_found = df_p.shape[0]
    ranks = int(df_p["clock_rank"].unique().max())
    if genes_found >= max_gene:
        print("Found interesting patient idx ",interesting_p, " with genes : ", genes_found)
        max_gene = genes_found
        df_int = df_p
    genes = {}
    for r in range(1,ranks+1):
        df_p_r =  df_p[df_p["clock_rank"]==r]
        #take Gene and Status 
        couple = df_p_r[['gene', 'mutatation_status']].drop_duplicates()
        genes_list = [f"{row[0]} ({row[1]})" for row in couple.values]        
        if len(genes_list)>0: #some ranks are skipped ????
            genes[r] = genes_list
    trees[p] = genes
    interesting_p += 1
    
'''for k,v in trees.items():
    for k2,v2 in v.items():   
        print(k,k2,v2)'''
    

# %%
from utils import plot_tree as pt


p = patients[0]

print(trees[p])

pt.print_tree(trees[p],"rank_same.gv")

# %% [markdown]
# ## Perchè lo stesso gene viene ripetuto? 
# Il gene è composto da tantissime basi (lettere); in questo dataset abbiamo l'elenco delle mutazioni (errori delle lettere), non dei geni. 
# 
# Quindi, se lo stesso gene compare più volte, vuol dire che in esso sono presenti più mutazioni puntiformi. Il gene si trova in uno specifico cromosoma e uno specifico segmento, l'unica cosa che cambia sono i dati di sequenziamento. 
# 
# DP/NV: riguardano come la sequenza del gene viene letta dal sequenziatore. 
# * DP: profondità, numero di volte in cui il sequenziatore ha letto quel segmento del dna; 
# * NV: numero di variants, numero di volte in cui leggendo quel punto, il sequenziatore ha trovato un variante (lettera sbagliata)
# 
# Quindi il tumore ha accumulato 5 mutazioni puntiformi diverse (M) in 5 punti diversi del gene STAG1. All'istante di clock_mean, il cromosoma che lo conteneva si è amplificato: loro non hanno trovato quando sono avvenute queste mutazioni, MA solo quando è avvenuta la duplicazione del segmento cromosomico(ovviamente tutte le 5 mutazioni avranno lo stesso clock_mean). 

# %% [markdown]
# 


