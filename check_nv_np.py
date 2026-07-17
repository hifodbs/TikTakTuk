import pyreadr
import pandas as pd

def check_multiplicity_data(df_filtrato):
    print("\n" + "="*50)
    print("--- 1. RICERCA NOME COLONNA MOLTEPLICITÀ ---")
    
    # Stampa tutte le colonne per aiutarti a trovare il nome corretto
    # print("Elenco di TUTTE le colonne nel dataset:")
    # print(df_filtrato.columns.tolist())
    
    col_molteplicita = 'mult_estimate' 
    
    print("\n" + "-"*50)
    if col_molteplicita not in df_filtrato.columns:
        print(f"ERRORE: La colonna '{col_molteplicita}' NON esiste nel dataset.")
        print("Guarda l'elenco delle colonne qui sopra, trova il nome giusto e modificalo nello script!")
        return
        
    print(f"[OK] Colonna '{col_molteplicita}' trovata con successo.")
    
    print("\n" + "-"*50)
    print("--- 2. ANALISI DEI DATI DI MOLTEPLICITÀ DOVE NV/DP SONO PRESENTI ---")
    
    # Filtriamo isolando SOLO le righe dove abbiamo i dati di lettura (NV e DP non sono NaN)
    mask_valid_reads = df_filtrato['NV'].notna() & df_filtrato['DP'].notna()
    df_valid = df_filtrato[mask_valid_reads]
    
    totale_valide = len(df_valid)
    print(f"Totale righe con dati NV e DP validi: {totale_valide}")
    
    if totale_valide > 0:
        # Tra queste righe valide, quante hanno il dato sulla molteplicità?
        mask_molt_presente = df_valid[col_molteplicita].notna()
        molt_presenti = mask_molt_presente.sum()
        molt_mancanti = totale_valide - molt_presenti
        
        print(f"-> Di queste, quante HANNO il dato '{col_molteplicita}'? {molt_presenti}")
        print(f"-> Quante NON lo hanno (sono NaN)? {molt_mancanti}")
        
        if molt_mancanti == 0:
            print("\nRISULTATO IDEALE: Hai il dato di molteplicità per tutte le mutazioni analizzabili! Puoi sostituire il VAF test manuale.")
        else:
            print(f"\nATTENZIONE: Ci sono {molt_mancanti} righe senza molteplicità.")
            
        print(f"\nDistribuzione dei valori di {col_molteplicita} (i primi 10 più frequenti):")
        print(df_valid[col_molteplicita].value_counts().head(10).to_string())

        # ==============================================================
        # NUOVA SEZIONE: RELAZIONE TRA CARIOTIPO E MOLTEPLICITÀ
        # ==============================================================
        print("\n" + "-"*50)
        print("--- 3. DISTRIBUZIONE MOLTEPLICITÀ PER CARIOTIPO ---")
        
        if 'karyotype' in df_valid.columns:
            # Prendiamo solo i dati che HANNO la molteplicità
            df_con_molt = df_valid[mask_molt_presente]
            
            # Raggruppiamo per Cariotipo e Molteplicità
            distribuzione = df_con_molt.groupby('karyotype')[col_molteplicita].value_counts().to_frame(name='conteggio').reset_index()
            # Ordiniamo i risultati per comodità di lettura
            distribuzione = distribuzione.sort_values(by=['karyotype', 'conteggio'], ascending=[True, False])
            
            print(f"\nConteggio delle stime di molteplicità (presenti) per ogni cariotipo:")
            print(distribuzione.to_string(index=False))
            
            # Analizziamo anche i mancanti!
            if molt_mancanti > 0:
                print(f"\n-> A quali cariotipi corrispondono le {molt_mancanti} righe con molteplicità MANCANTE (NaN)?")
                df_senza_molt = df_valid[~mask_molt_presente]
                print(df_senza_molt['karyotype'].value_counts().to_string())
        else:
            print("\nImpossibile fare l'analisi: la colonna 'karyotype' non è presente.")


    print("\n" + "="*50 + "\n")


def run_global_check():
    filepath = '06_Cb_BTM_table.rds'
    print(f"Caricamento intero dataset da {filepath}...")
    
    # Carica l'intero dataset
    result = pyreadr.read_r(filepath)
    df = result[None]
    
    # Hai rimosso il filtro WT, quindi processiamo tutto il dataframe
    check_multiplicity_data(df)

if __name__ == "__main__":
    run_global_check()