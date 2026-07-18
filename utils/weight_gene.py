import numpy as np
import pandas as pd


def weight_genes_bypass(df_filtrato):
    df = df_filtrato.copy()
    df = df[["gene", "sample_id", "clock_rank"]]

    # Keep only unique (gene, sample, rank) records
    df = df.drop_duplicates()

    # Count the number of unique rows per gene and sample
    df_counts = (
        df.groupby(["gene", "sample_id"]).size().reset_index(name="count")
    )

    return df_counts


def weight_genes(df_filtrato, top_n=20, bypass=False):
    if bypass:
        return weight_genes_bypass(df_filtrato)
    df = df_filtrato.copy()

    # 1. Mappatura tipologia mutazione (WT neutro a 1.0 per chi non ha WGD)
    status_weights = {'CI_M': 2.0, 'M': 2.0, 'CNA_driver': 4.0, 'WT': 1.0}
    df['w_status'] = df['mutatation_status'].map(status_weights).fillna(1.0)

    # 2. Moltiplicatore PCAWG
    df['b_pcawg'] = np.where(df['mutation_call'].notna(), 2.0, 1.0)

    # 3. NUOVA LOGICA CARIOTIPO / MOLTEPLICITÀ
    # Se mult_estimate è presente (non è NaN), usiamo 1.0 + (molteplicità * 0.2), cariotipo: 2:0, 2:2, 2:1
    # Se mult_estimate è NaN, usiamo 1.0: WT, CNA_driver
    df['k_severity'] = np.where(
        df['mult_estimate'].notna(), 
        1.0 + (df['mult_estimate'] * 0.2), 
        1.0
    )

    # 4. Fattore WGD globale
    df['wgd_factor'] = np.where(df['is_WGD'] == True, 0.8, 1.0)

    # 5. Calcolo peso della riga
    df['row_weight'] = (df['w_status'] * df['b_pcawg'] * df['k_severity'] * df['wgd_factor'])
    
    # Ordiniamo per isolare la mutazione più grave
    df = df.sort_values(by=['clock_rank', 'gene', 'row_weight'], ascending=[True, True, False])

    # 6. Filtro pazienti multipli (1 paziente = 1 voto col peso peggiore)
    df_pazienti_unici = df.groupby(['sample_id', 'clock_rank', 'gene']).first().reset_index()

    # 7. Aggregazione per RANK
    node_weights = df_pazienti_unici.groupby(['clock_rank', 'gene']).agg(
        gene_weight=('row_weight', 'sum')
    ).reset_index()
    
    # 8. ESTRAZIONE TOP 20 PER RANK
    top_genes_df = node_weights.sort_values(
        by=['clock_rank', 'gene_weight'], ascending=[True, False]
    ).groupby('clock_rank').head(top_n)
    
    # Creiamo un dizionario che associa a ogni Rank i suoi 20 geni: {1: {'TP53', ...}, 2: {'PIK3CA', ...}}
    valid_rank_genes = top_genes_df.groupby('clock_rank')['gene'].apply(set).to_dict()
    
    # Pesi globali per il clustering (facciamo la media sui Top estratti)
    final_weights = top_genes_df.groupby('gene')['gene_weight'].mean().to_dict()
    
    return final_weights, valid_rank_genes