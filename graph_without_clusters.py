import pyreadr
import pandas as pd
import collections
import math
from utils import plot_tree as pt

# ==========================================
# 1. PARAMETRI E CARICAMENTO DATI
# ==========================================
ttype_selezionato = "OV"
min_pazienti_condivisi = 0  # Aumenta questo valore per filtrare archi rari

print(f"Calcolo delle co-occorrenze per {ttype_selezionato} in corso...")
result = pyreadr.read_r('06_Cb_BTM_table.rds')
df = result[None]

# Filtraggio
df_filtrato = df[(df["ttype"] == ttype_selezionato) & (df['mutatation_status'] != "WT")].copy()
df_filtrato = df_filtrato.dropna(subset=["clock_rank", "gene"])

# ==========================================
# 2. ESTRAZIONE TRAIETTORIE DEI PAZIENTI
# ==========================================
patient_genes = collections.defaultdict(lambda: collections.defaultdict(set))

# Per ogni riga valida, aggiungiamo il gene al corrispondente rank del paziente
for _, row in df_filtrato.iterrows():
    p_id = row['sample_id']
    rank = int(row['clock_rank'])
    gene = row['gene']
    patient_genes[p_id][rank].add(gene) # set() evita automaticamente che lo stesso gene venga contato 2 volte

# ==========================================
# 3. CALCOLO CO-OCCORRENZE E FREQUENZE
# ==========================================
co_occurrences = collections.defaultdict(int)
gene_freq_source = collections.defaultdict(int)
gene_freq_target = collections.defaultdict(int)
total_transitions = 0

for p_id, ranks_dict in patient_genes.items():
    active_ranks = sorted(ranks_dict.keys())
    
    # Se il tumore ha avuto una sola ondata, non c'è co-occorrenza temporale da tracciare
    if len(active_ranks) < 2:
        continue
        
    for i in range(len(active_ranks) - 1):
        r_source = active_ranks[i]
        r_target = active_ranks[i + 1]
        
        genes_source = ranks_dict[r_source]
        genes_target = ranks_dict[r_target]
        
        # Frequenze marginali
        for g_s in genes_source:
            gene_freq_source[g_s] += 1
        for g_t in genes_target:
            gene_freq_target[g_t] += 1
            
        # Conteggio congiunto: quante volte A è seguito da B
        for g_s in genes_source:
            for g_t in genes_target:
                co_occurrences[(r_source, g_s, r_target, g_t)] += 1
                total_transitions += 1


# ==========================================
# 4. CALCOLO PMI E CREAZIONE STRUTTURA GRAFO
# ==========================================
valid_nodes = collections.defaultdict(set)

# INSERIMENTO DI TUTTI I GENI NEL GRAFO (Anche quelli isolati)
for p_id, ranks_dict in patient_genes.items():
    for rank, genes in ranks_dict.items():
        for g in genes:
            valid_nodes[rank].add(g)

edges = []

if total_transitions > 0:
    for (r_s, g_s, r_t, g_t), count in co_occurrences.items():
        p_ab = count / total_transitions
        p_a = gene_freq_source[g_s] / total_transitions
        p_b = gene_freq_target[g_t] / total_transitions
        
        # PMI Formula
        pmi = math.log2(p_ab / (p_a * p_b))
        
        # Filtro: PMI positiva e minimo di pazienti
        if pmi > 0 :
            edges.append({
                'source': f"{g_s}_R{r_s}",
                'target': f"{g_t}_R{r_t}",
                'weight': count,
                'pmi': pmi
            })

# Conversione dei set in liste per comodità
for rank in valid_nodes:
    valid_nodes[rank] = list(valid_nodes[rank])

# ==========================================
# 5. GENERAZIONE GRAFO
# ==========================================
file_name = f"gene_pmi_tree_{ttype_selezionato}.gv"
print(f"Calcolo terminato. Generazione albero: {file_name}")

pt.print_tree_genes_pmi(valid_nodes, edges, file_name)