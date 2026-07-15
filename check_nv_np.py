import pyreadr
import pandas as pd

def check_nv_dp_data_updated(df_filtrato):
    print("\n" + "="*50)
    print("--- 1. ANALISI DEI VALORI MANCANTI (NaN) ---")
    
    # Filtriamo solo le righe dove NV è NaN
    nan_mask = df_filtrato['NV'].isna()
    df_nan = df_filtrato[nan_mask]
    
    print(f"Totale righe con NV/DP mancanti: {len(df_nan)}")
    
    if len(df_nan) > 0:
        print("\n-> Distribuzione 'mutatation_status' nei mancanti:")
        print(df_nan['mutatation_status'].value_counts().to_string())
            
        if 'karyotype' in df_nan.columns:
            # Nel formato Major:Minor, la LOH si identifica con l'allele minore a 0 (termina con ':0')
            loh_mask = df_nan['karyotype'].astype(str).str.endswith(':0')
            loh_count = loh_mask.sum()
            print(f"\n-> Quanti dei {len(df_nan)} mancanti sono in stato di LOH (terminano con ':0')? {loh_count}")
            
            print("\n-> Dettaglio esatto dei cariotipi per queste righe mancanti:")
            print(df_nan['karyotype'].value_counts().to_string())

    print("\n" + "-"*50)
    print("--- 2. ANALISI VALORI MULTIPLI PER SINGOLO PAZIENTE ---")
    
    if 'segment_id' in df_filtrato.columns and 'sample_id' in df_filtrato.columns:
        # Raggruppamento inserendo il sample_id (Paziente singolo)
        stats_paz = df_filtrato.groupby(['sample_id', 'gene', 'segment_id']).agg(
            nv_unici=('NV', 'nunique'),
            dp_unici=('DP', 'nunique'),
            numero_occorrenze=('NV', 'size')
        ).reset_index()
        
        # Filtriamo per trovare anomalie all'interno dello stesso paziente
        variazioni_paz = stats_paz[(stats_paz['nv_unici'] > 1) | (stats_paz['dp_unici'] > 1)]
        
        print(f"Combinazioni uniche (Paziente + Gene + Segmento): {len(stats_paz)}")
        print(f"Casi in cui un SINGOLO PAZIENTE ha valori NV/DP multipli: {len(variazioni_paz)}")
        
        if len(variazioni_paz) > 0:
            print("\n-> Esempio (i primi 5 casi anomali per singolo paziente):")
            print(variazioni_paz.head().to_string())
        else:
            print("\n-> PERFETTO! A livello di singolo paziente, NV e DP sono sempre univoci per segmento.")
    else:
        print("\nColonne 'sample_id' o 'segment_id' non trovate.")
        
    print("="*50 + "\n")

def run_global_check():
    filepath = '06_Cb_BTM_table.rds'
    print(f"Caricamento intero dataset da {filepath}...")
    
    # Carica l'intero dataset
    result = pyreadr.read_r(filepath)
    df = result[None]
    
    # Filtro: rimuoviamo i WT per analizzare solo le mutazioni
    df_mutated = df[df['mutatation_status'] != "WT"].copy()
    
    # Lancia il controllo
    check_nv_dp_data_updated(df_mutated)

if __name__ == "__main__":
    run_global_check()