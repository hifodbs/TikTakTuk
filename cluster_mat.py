import pyreadr
import numpy as np
import pandas as pd
import math
from utils import plot_tree as pt
from utils.weight_gene import weight_genes
import collections
from utils.cluster import clusterize_patients
import utils.plot_both
import networkx as nx
import matplotlib.pyplot as plt
import os
#import utils.plot_cancer_progession

#-----------------------------------
#Caricamento e pulizia dei dati 
#----------------------------------

def load_and_filter_data(filepath, ttype_selezionato):
    """Carica i dati e filtra per tipo di tumore e geni validi."""
    
    print(f"Caricamento dataset da {filepath}...")
    result = pyreadr.read_r(filepath)
    df = result[None]

    if ttype_selezionato != None:
        df_filtrato = df[df["ttype"] == ttype_selezionato].copy()  #selezione del tumore 
    else:
        df_filtrato = df.copy()  #selezione del tumore 
    #esclusione wgd-wt
    is_wgd_wt = (df_filtrato['mutatation_status'] == 'WT') & (df_filtrato['is_WGD'] == True)
    df_filtrato = df_filtrato[~is_wgd_wt]
    df_filtrato = df_filtrato.dropna(subset=["clock_rank", "gene"])
    df_filtrato.reset_index(drop=True, inplace=True)
    
    print(f"Dataset filtrato per {ttype_selezionato}. Righe: {len(df_filtrato)}")
    print("number of patients",len(df["sample_id"].unique()))
    return df_filtrato

def print_patient_ranks(df):
    """Stampa tutti i rank associati a ciascun paziente."""
    patient_ranks = collections.defaultdict(list)

    for _, row in df.iterrows():
        patient_ranks[row["sample_id"]].append(int(row["clock_rank"]))

    for sample_id in sorted(patient_ranks):
        unique_ranks = sorted(set(patient_ranks[sample_id]))
        print(f"Patient {sample_id}: ranks {unique_ranks}")


#-----------------------------------------------------
#Assegnamento del peso ad ogni gene
#e mappatura ID di ogni paziente ad un numero intero
#----------------------------------------------------

def build_boolean_matrix(df_filtrato, top_n=20):
    """Calcola i pesi, filtra per Rank e costruisce l'array booleano omogeneo."""
    import collections
    input_dict = collections.defaultdict(dict)
    
    # 1. Otteniamo i pesi e la mappa dei Top 20 per Rank
    weights, valid_rank_genes = weight_genes(df_filtrato, top_n)
    
    # all_genes sarà l'unione di tutti i Top 20 (serve per far funzionare Jaccard nel clustering)
    all_genes = list(weights.keys())
    all_weights = list(weights.values())
    
    print(f"Costruzione Matrice: Trovati {len(all_genes)} geni unici unendo i Top {top_n} di ogni Rank.")
    
    df_subset = df_filtrato[["sample_id", "gene", "clock_rank"]]
    
    for sample_id in df_subset["sample_id"].unique():
        df_sample = df_subset[df_subset["sample_id"] == sample_id]
        
        for rank in df_sample["clock_rank"].unique():
            df_s_rank = df_sample[df_sample["clock_rank"] == rank]
            patient_genes = df_s_rank["gene"].to_numpy()
            
            # FILTRO CHIAVE: Isoliamo solo i geni ammessi per QUESTO specifico Rank
            valid_for_this_rank = valid_rank_genes.get(rank, set())
            patient_genes_filtered = [g for g in patient_genes if g in valid_for_this_rank]
            
            # Array booleano
            result = np.isin(all_genes, patient_genes_filtered)
            input_dict[rank][sample_id] = result
            
    return input_dict, all_genes, all_weights

def get_patient_mapping(input_dict):
    """Crea un dizionario che associa ogni UUID a un numero intero."""
    all_p_ids = set()
    for rank in input_dict:
        for p_id in input_dict[rank].keys():
            all_p_ids.add(p_id)
            
    # Ordiniamo gli ID per garantire che la numerazione sia sempre la stessa a ogni avvio
    mapping = {p_id: str(i + 1) for i, p_id in enumerate(sorted(list(all_p_ids)))}
    return mapping

#-----------------------------------------------
#clustering
#----------------------------------------------
def run_clustering(input_dict, all_weights, distance_metric="jaccard", linkage="average", best_k_metric="silhouette"):
    """Esegue il clustering gerarchico sui pazienti."""
    distance_matrixes, labelsels = clusterize_patients(
        input_dict, distance_metric, linkage, best_k_metric, all_weights
    )
    return distance_matrixes, labelsels

# linkage = "average"           #metodo per unire i cluster, calcolando la dist media tra i membri dei due clusters
# distance_metric = "jaccard"   #measure the similarity 
# best_k_metric = "silhouette"  #divide i pazienti in un numero diverso di clusters e in base a Sil. trova il k migliore


#-----------------------------------------------
#associazione paziente ai clusters
#----------------------------------------------

def extract_patients_clusters(input_dict, labelsels, patient_mapping):
    """Associa i pazienti (ID accorciato per leggibilità) ai cluster assegnati."""
    patient_clusters = collections.defaultdict(lambda: collections.defaultdict(list))
    
    for rank in input_dict.keys():
        patient_ids = list(input_dict[rank].keys())
        labels = labelsels[rank]
        
        for idx, p_id in enumerate(patient_ids):
            c_label = labels[idx]
            # Id del paziente in numero intero
            num_id = patient_mapping[p_id]
            patient_clusters[rank][c_label].append(num_id)
            
    return patient_clusters

#--------------------------------------------------
#associaizone geni ai cluster
#--------------------------------------------------

def extract_genes_clusters(input_dict, labelsels, all_genes, patient_mapping):
    """Associa i geni ai cluster salvando l'elenco dei numeri dei pazienti che li possiedono."""
    # Struttura più sicura per evitare KeyError
    gene_clusters = collections.defaultdict(lambda: collections.defaultdict(dict))
    
    for rank in input_dict.keys():
        patient_ids = list(input_dict[rank].keys())
        labels = labelsels[rank]
        
        for idx, p_id in enumerate(patient_ids):
            c_label = labels[idx]
            num_id = patient_mapping[p_id]
            bool_arr = input_dict[rank][p_id]
            
            for i, is_mutated in enumerate(bool_arr):
                if is_mutated:
                    gene = all_genes[i]
                    
                    # Se il gene non c'è ancora in questo cluster, inizializza una lista vuota
                    if gene not in gene_clusters[rank][c_label]:
                        gene_clusters[rank][c_label][gene] = []
                        
                    # Aggiunge il numero del paziente alla lista del gene
                    gene_clusters[rank][c_label][gene].append(num_id)
                        
    return gene_clusters

#-----------------------------------------------
#plots dei clusters
#----------------------------------------------

#funzioni separate per scatter plots e dendogramma
#utils.plot_patient_clustered.plot(input, labelsels,distance_matrixes)

#import utils.plot_patients_dendograms

#utils.plot_patients_dendograms.plot(input, labelsels,distance_matrixes)

#funzione per scatter plots e dendogramma insieme 
#utils.plot_both.plot(input, labelsels,distance_matrixes)



#-------------------------------------------------------------------
#find the important edges
#------------------------------------------------------------------

def extract_patient_trajectories(input_dict, labelsels, patient_mapping):
    """Estrae la mappa dei salti di ogni paziente: {num_paziente: {rank: cluster_id}}"""
    trajectories = collections.defaultdict(dict)
    
    for rank in input_dict.keys():
        patient_ids = list(input_dict[rank].keys())
        labels = labelsels[rank]
        
        for idx, p_id in enumerate(patient_ids):
            num_id = patient_mapping[p_id]
            c_label = labels[idx]
            trajectories[num_id][rank] = c_label
            
    return trajectories

def calculate_cluster_pmi(trajectories, min_pazienti=1):
    """Calcola la Mutual Information (PMI) per validare le transizioni tra i cluster."""
    import math
    cluster_transitions = collections.defaultdict(int)
    cluster_freq_source = collections.defaultdict(int)
    cluster_freq_target = collections.defaultdict(int)
    total_transitions = 0

    # 1. Conteggio frequenze assolute
    for p_num, path in trajectories.items():
        active_ranks = sorted(path.keys())
        if len(active_ranks) < 2:
            continue
            
        for i in range(len(active_ranks) - 1):
            r_s = active_ranks[i]
            c_s = path[r_s]
            r_t = active_ranks[i + 1]
            c_t = path[r_t]
            
            source_id = f"R{r_s}_C{c_s}"
            target_id = f"R{r_t}_C{c_t}"
            
            cluster_freq_source[source_id] += 1
            cluster_freq_target[target_id] += 1
            cluster_transitions[(r_s, c_s, r_t, c_t)] += 1
            total_transitions += 1

    # 2. Calcolo probabilità e filtraggio PMI
    pmi_edges = []
    if total_transitions > 0:
        for (r_s, c_s, r_t, c_t), count in cluster_transitions.items():
            source_id = f"R{r_s}_C{c_s}"
            target_id = f"R{r_t}_C{c_t}"
            
            p_ab = count / total_transitions
            p_a = cluster_freq_source[source_id] / total_transitions
            p_b = cluster_freq_target[target_id] / total_transitions
            
            pmi = math.log2(p_ab / (p_a * p_b))
            
            # Manteniamo la freccia solo se la PMI è positiva (dipendenza reale)
            if pmi > 0 and count >= min_pazienti:
                pmi_edges.append({
                    'source': (r_s, c_s),
                    'target': (r_t, c_t),
                    'weight': count,
                    'pmi': pmi
                })
                
    return pmi_edges

# GLOBAL 

def calculate_cluster_pmi_global(input_dict, labelsels, patient_mapping):
    import math
    import collections
    
    # 1. Mappiamo correttamente: Rank -> {Num_Paziente: ID_Cluster}
    # Questo risolve il problema dello slittamento degli indici
    rank_patient_cluster = collections.defaultdict(dict)
    
    for rank in input_dict.keys():
        patient_ids = list(input_dict[rank].keys())
        labels = labelsels[rank]
        for idx, p_id in enumerate(patient_ids):
            num_id = patient_mapping[p_id]
            rank_patient_cluster[rank][num_id] = labels[idx]

    all_ranks = sorted(rank_patient_cluster.keys())
    pmi_edges = []
    
    # 2. Troviamo il numero totale di pazienti univoci (N globale) 
    # Serve per calcolare probabilità veritiere (P(A), P(B))
    all_patients = set()
    for rank in rank_patient_cluster:
        all_patients.update(rank_patient_cluster[rank].keys())
    N = len(all_patients)
    
    # 3. Confrontiamo tutte le combinazioni possibili per Rank consecutivi
    for i in range(len(all_ranks) - 1):
        r1 = all_ranks[i]
        r2 = all_ranks[i + 1]
        
        # Set dei cluster univoci presenti in ciascun Rank
        clusters_r1 = set(rank_patient_cluster[r1].values())
        clusters_r2 = set(rank_patient_cluster[r2].values())
        
        for c1 in clusters_r1:
            for c2 in clusters_r2:
                # Estraiamo i set di pazienti che appartengono a c1(Rank 1) e c2(Rank 2)
                p_in_c1 = {p for p, c in rank_patient_cluster[r1].items() if c == c1}
                p_in_c2 = {p for p, c in rank_patient_cluster[r2].items() if c == c2}
                
                # Quanti pazienti hanno fatto questo esatto percorso?
                co_occ = len(p_in_c1.intersection(p_in_c2))
                
                # Calcoliamo la PMI solo se c'è almeno 1 paziente in comune per evitare log(0)
                if co_occ > 0:
                    p_ab = co_occ / N          # Probabilità congiunta
                    p_a = len(p_in_c1) / N     # Probabilità del cluster 1
                    p_b = len(p_in_c2) / N     # Probabilità del cluster 2
                    
                    pmi = math.log2(p_ab / (p_a * p_b))
                    
                    # Filtriamo tenendo solo correlazioni positive (maggiore di 0)
                    if pmi > 0: 
                        pmi_edges.append({
                            'source': (r1, c1),
                            'target': (r2, c2),
                            'pmi': pmi
                        })
                        
    return pmi_edges

def create_gene_connections(sorted_pmi_edges, gene_clusters):
    genes = []
    for e in sorted_pmi_edges:
        for start_gene in gene_clusters[e["source"][0]][e["source"][1]].keys():
            for end_gene in gene_clusters[e["target"][0]][e["target"][1]].keys():
                genes.append((start_gene,end_gene,e["pmi"]))
    return  sorted(genes, key=lambda x: x[2],reverse=True)

def create_tree(gene_connections):
    graph = nx.DiGraph()
    
    for conn in gene_connections:
        if graph.has_node(conn[0]) and graph.has_node(conn[1]) and nx.has_path(graph,conn[1],conn[0]):
            continue
        graph.add_edge(conn[0],conn[1], weight=conn[2])
        
    return graph
    
def main():
    # Parametri
    file_rds = '06_Cb_BTM_table.rds'
    tumore = "ESAD" #None
    
    # 1. Caricamento e Preparazione
    df = load_and_filter_data(file_rds, tumore)
    input_dict, all_genes, all_weights = build_boolean_matrix(df, top_n=20)
    patient_mapping = get_patient_mapping(input_dict)
    
    # 2. Esecuzione Clustering
    print("Avvio clustering dei pazienti...")
    distance_matrixes, labelsels = run_clustering(input_dict, all_weights)
    
    # 3. Estrazione dei Dizionari Dati
    patient_clusters = extract_patients_clusters(input_dict, labelsels, patient_mapping)
    gene_clusters = extract_genes_clusters(input_dict, labelsels, all_genes, patient_mapping)
    
    # 4. Calcolo Traiettorie e PMI per l'evoluzione
    print("Calcolo Mutual Information (PMI) tra cluster...")
    trajectories = extract_patient_trajectories(input_dict, labelsels, patient_mapping)
    pmi_edges = calculate_cluster_pmi(trajectories, min_pazienti=1)
    
    # 5. Generazione Archi Semplici (tutto Rank) per visualizzazione di controllo
    rank_edges = pt.generate_rank_edges(input_dict)
    
#    --- STAMPE GRAFICHE ---
    # A. Albero semplificato (Pazienti)
    print("Generazione albero Pazienti semplice...")
    pt.plot_patients_clusters_simple(patient_clusters, rank_edges, f"tree_patients_simple_{tumore}.gv")
    
#    B. Albero semplificato (Geni)
    print("Generazione albero Geni semplice...")
    pt.plot_genes_clusters_simple(gene_clusters, rank_edges, f"tree_genes_simple_{tumore}.gv")
    
    # C. Albero evolutivo reale basato sulle frequenze (PMI)
    print("Generazione albero evolutivo PMI...")
    pt.plot_genes_clusters_pmi(gene_clusters, pmi_edges, f"tree_genes_PMI_{tumore}.gv")
    
    # sorted_pmi_edges = sorted(pmi_edges, key=lambda x: x['pmi'])
    
    # gene_connections = create_gene_connections(sorted_pmi_edges, gene_clusters)
    # # Keeps only items with pmi >= 2.5
    # filtered_gene_connection = [d for d in gene_connections if d[2] >= 2.5]
    
    
    # gad = create_tree(filtered_gene_connection)
    # app = utils.plot_cancer_progession.create_interactive_dag_app(gad)
    # app.run(debug=True)


if __name__ == "__main__":
    main()