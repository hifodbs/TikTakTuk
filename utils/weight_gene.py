import numpy as np
import pandas as pd

import numpy as np
import pandas as pd

def check_if_amplified_or_loh(row):
    """
    Classifica la singola riga coprendo i tre scenari:
    1. LOH (restituisce True)
    2. NaN senza LOH (restituisce False)
    3. Dati presenti (calcola la distanza VAF, True se amplificato)
    """
    karyo = str(row['karyotype'])
    try:
        maj, min_allele = map(int, karyo.split(':'))
    except ValueError:
        return False
        
    if min_allele > maj:
        maj, min_allele = min_allele, maj

    cn = maj + min_allele

    # ==========================================
    # CASO 1: LOH (Loss of Heterozygosity) o Amplificazione Massiccia
    # ==========================================
    if min_allele == 0 or cn >= 5:
        return True 
        
    # ==========================================
    # CASO 2: Valori NaN (ma senza LOH)
    # ==========================================
    if pd.isna(row['NV']) or pd.isna(row['DP']) or row['DP'] == 0:
        return False

    # ==========================================
    # CASO 3: Dati presenti, calcolo della VAF
    # ==========================================
    if maj == min_allele or cn == 0:
        return False  # Cariotipo bilanciato
        
    vaf = row['NV'] / row['DP']
    e_maj = maj / cn
    e_min = min_allele / cn
    
    # Test della distanza
    dist_maj = abs(vaf - e_maj)
    dist_min = abs(vaf - e_min)
    
    # True se la mutazione appartiene all'allele amplificato (Major)
    return dist_maj < dist_min


def weight_genes(df_filtrato, top_n=20):
    df = df_filtrato.copy()

    # 1. Mappatura tipologia mutazione (WT neutro a 1.0 per chi non ha WGD)
    status_weights = {'CI_M': 2.0, 'M': 2.0, 'CNA_driver': 4.0, 'WT': 1.0}
    df['w_status'] = df['mutatation_status'].map(status_weights).fillna(1.0)

    # 2. Moltiplicatore PCAWG
    df['b_pcawg'] = np.where(df['mutation_call'].notna(), 2.0, 1.0)

    # 3. Valutazione Cariotipo (Richiede la tua funzione check_if_amplified_or_loh)
    df['is_severe_karyo'] = df.apply(check_if_amplified_or_loh, axis=1)
    df['k_severity'] = np.where(df['is_severe_karyo'], 1.5, 1.0)

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