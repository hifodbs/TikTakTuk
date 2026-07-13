import numpy as np
import pandas as pd


def weight_genes(df_filtrato):
    # funzioni per la mappatura dei pesi

    # Mappatura tipologia mutazione per distaccare i driver principali
    status_weights = {'CI_M': 2.0, 'M': 2.0, 'CNA_driver': 4.0}
    df_filtrato['w_status'] = df_filtrato['mutatation_status'].map(status_weights).fillna(1.0)

    # Moltiplicatore PCAWG (2.0 se il catalogo conferma il driver, 1.0 altrimenti)
    df_filtrato['b_pcawg'] = np.where(df_filtrato['mutation_call'].notna(), 2.0, 1.0)

    # Moltiplicatore Cariotipo (LOH o instabilità numerica)
    def analizza_cariotipo(k_str):
        if pd.isna(k_str): 
            return 1.0
        try:
            parts = str(k_str).split(':')
            major, minor = int(parts[0]), int(parts[1])
            if minor == 0: 
                return 1.5  # Incremento per Loss of Heterozygosity
            if (major + minor) >= 5: 
                return 1.5  # Incremento per amplificazione massiccia
        except:
            pass
        return 1.0

    df_filtrato['k_severity'] = df_filtrato['karyotype'].apply(analizza_cariotipo)

    # Fattore WGD globale del paziente
    df_filtrato['wgd_factor'] = np.where(df_filtrato['is_WGD'] == True, 0.8, 1.0)

    # Calcolo del peso finale della singola riga (istanza mutazionale)
    df_filtrato['gene_weight'] = (df_filtrato['w_status'] * df_filtrato['b_pcawg'] * df_filtrato['k_severity'] * df_filtrato['wgd_factor'])
    
    df_filtrato = df_filtrato.sort_values(by=['clock_rank', 'gene', 'w_status'], ascending=[True, True, False])

    # Aggregation
    # Somma automaticamente le mutazioni multiple dello stesso gene nello stesso rank e tra diversi pazienti
    node_weights = df_filtrato.groupby(['clock_rank', 'gene']).agg(
        gene_weight=('gene_weight', 'sum'),                  # Somma i pesi per il punteggio totale
        mut_status=('mutatation_status', 'first'),           # Conserva la tipologia più grave (grazie all'ordinamento di prima)
        n_pazienti=('sample_id', 'nunique'),                 # Conta in quanti pazienti DIVERSI compare
        n_occorrenze=('gene', 'count')                       # Conta il numero totale di mutazioni subite in assoluto
    ).reset_index()
    
    return pd.Series(node_weights.gene_weight.values,index=node_weights.gene).to_dict()