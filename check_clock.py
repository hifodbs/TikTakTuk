import pyreadr
import pandas as pd

def check_full_distribution():
    filepath = '06_Cb_BTM_table.rds'
    print(f"Caricamento dataset da {filepath}...")
    
    result = pyreadr.read_r(filepath)
    df = result[None]
    
    col_molt = 'mult_estimate'
    col_status = 'mutatation_status' 
    col_clock = 'clock_rank'
    col_wgd = 'is_WGD'
    
    # Filtriamo tenendo solo le righe con un CLOCK_RANK valido
    df_with_clock = df[df[col_clock].notna()].copy()
    
    # Creiamo i due gruppi
    df_con_molt = df_with_clock[df_with_clock[col_molt].notna()].copy()
    df_senza_molt = df_with_clock[df_with_clock[col_molt].isna()].copy()
    
    print("\n" + "="*80)
    print(f"DISTRIBUZIONE TOTALE: {len(df_with_clock)} righe con clock_rank valido")
    print("="*80)
    
    print(f"\n1. GRUPPO CON MOLTEPLICITÀ ({len(df_con_molt)} righe):")
    tabella_con = pd.crosstab(df_con_molt[col_status], df_con_molt[col_wgd], margins=True, margins_name="Totale")
    print(tabella_con.to_string())
    
    print("\n" + "-"*50)
    
    print(f"\n2. GRUPPO SENZA MOLTEPLICITÀ ({len(df_senza_molt)} righe):")
    tabella_senza = pd.crosstab(df_senza_molt[col_status], df_senza_molt[col_wgd], margins=True, margins_name="Totale")
    print(tabella_senza.to_string())
    
    print("\n" + "="*80 + "\n")

if __name__ == "__main__":
    check_full_distribution()