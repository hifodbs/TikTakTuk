import pyreadr
import pandas as pd
import numpy as np
import math
from collections import defaultdict, Counter
import os
from utils import plot_tree as pt
import webbrowser

# Caricamento dei dati dal file RDS
print("Caricamento dataset in corso...")
result = pyreadr.read_r('06_Cb_BTM_table.rds')
df = result[None]

# Filtraggio
ttype_selezionato = "ESAD"  #selezione del tumore
df_filtrato = df[df["ttype"] == ttype_selezionato].copy()

# Eliminazione geni sani (WT) e dati mancanti critici
df_filtrato = df_filtrato[df_filtrato['mutatation_status'] != "WT"]
df_filtrato = df_filtrato.dropna(subset=["clock_rank", "gene"])

print(f"Dataset filtrato per {ttype_selezionato}. Righe totali da analizzare: {len(df_filtrato)}")

df = df_filtrato[["sample_id","gene","clock_rank"]]

df.reset_index(drop=True, inplace=True)



print("number of patients",len(df["sample_id"].unique()))

import utils.weight_gene
import collections

input = collections.defaultdict(dict)

weights = utils.weight_gene.weight_genes(df_filtrato)
all_genes = list(weights.keys())
all_weights = list(weights.values())

for sample_id in df["sample_id"].unique():
    df_sample = df[df["sample_id"]==sample_id]
    #print(sample_id," con ranks ", df_sample["clock_rank"].unique())
    for rank in df_sample["clock_rank"].unique():
        df_s_rank = df_sample[df_sample["clock_rank"] == rank]
        genes = df_s_rank["gene"].to_numpy()
        result = np.isin(all_genes,genes)
        input[rank][sample_id]=result

import utils.plot_patient_clustered
from utils.cluster import clusterize_patients

linkage = "average"
distance_metric = "jaccard"
best_k_metric = "silhouette"


distance_matrixes,labelsels = clusterize_patients(input, distance_metric, linkage, best_k_metric, all_weights)


#utils.plot_patient_clustered.plot(input, labelsels,distance_matrixes)


#import utils.plot_patients_dendograms

#utils.plot_patients_dendograms.plot(input, labelsels,distance_matrixes)

import utils.plot_both

utils.plot_both.plot(input, labelsels,distance_matrixes)
