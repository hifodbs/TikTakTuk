# %%
import pyreadr
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# %% [markdown]
# Il dataset mappa le alterazioni strutturali del numero di copie (CNA):
# - amplificazioni genomiche: trisomie (2:1) o tetrasomie 2:2 - aggiunta di un allele nel cromosoma
# - perdita di eterozigosi (LOH) - perdita di un allele nel cromosoma
# 
# TickTack è un modello statistico che prende i dati di sequenziamento dei tumori e:
# - calcola il tempo molecolare - dei segmenti cromosomici(tau): inserisce ogni amplificazione su una linea temporale da 0 a 1 (nascita/momento del prelievo del campione)
# - raggruppa gli eventi simultanei (Clustering dei Co-occorrenti): segmenti di cromosomi che si sono amplificati o diminuiti, nello stesso momenti; il numero dei macro-eventi temporali in un paziente è best_k.

# %%
# Load the rds file
result = pyreadr.read_r('06_Cb_BTM_table.rds')

# pyreadr returns a dictionary; the data is usually under the 'None' key 
# because rds files contain a single object without a specific name.
df = result[None]

# Show all columns

print("-"*50+"Columns name"+50*"-")
print(df.columns)
print("-"*50+"Columns name"+50*"-")
df


# %% [markdown]
# Il dataset è composto da 89696 righe che corrispondono a singoli geni. Questi geni sono stati inseriti nella tabella perché si trovavano a bordo di segmenti che hanno subìto un'alterazione strutturale (un'amplificazione o una perdita di eterozigosi copia-neutra), la maggior parte sono sani, solo alcuni hanno una mutazione.

# %%
# Show values of category columns
print("-"*50+"Labels values"+50*"-")
for col in df.columns:
    if(len(df[col].unique())<100 ):
        print(col + ": " + str(df[col].unique()))
print("-"*50+"Labels values"+50*"-")
print()

# %% [markdown]
# Explanation: 
# 1. Identificatori e geni: 
#     1.  segment_id: La regione cromosomica (segmento) in cui si trova la mutazione (i tumori perdono o duplicano interi pezzi di cromosomi):
#     - ex.: chr1:33554784:48120521:2:2:1
#         1. chr1: cromosoma 1
#         2. 33554784: Start position- punto di inizio del segmento 
#         3. 48120521: End position- punto finale del segmento 
#         4. 2:2:1  : karyotype (4 coppie totali, 2:2) e stato clonale (1)
#         
#     Se si guarda l'albero genealogico, si possono distinguere: 
#     - tronco dell'albero: stato clonale (ID 1)= qualsiasi mutazione che avviene in questa fase, verrà ereditata da tutti i rami e da tutte le foglie-   sarà presente nel 100% delle cellule tumorali
#     - rami dell'albero: stato subclonale (ID 2,3, ect..) - il danno genetico appartiene solo a una parte del tumore
#     2.  gene: Il gene in cui è stata identificata la mutazione
#     3. sample_id: identificativo del paziente e quindi tessuto tumorale
#     - ex.: 0009b464-b376-4fbc-8a56-da538269a02f: UUID (Universally Unique Identifier)
# 
# 2. Dati di sequenziamento:
#     1. karyotype: Il cariotipo locale del segmento genomico espresso come "Allele Maggiore : Allele Minore". Ad esempio, "2:1" indica due copie materne e una paterna, "2:0" indica una perdita di eterozigosi (LOH), e "2:2" un assetto bilanciato.
#     2. NV (Variant Nucleotide Reads): Il numero di letture (reads) di sequenziamento che supportano la mutazione.
#     3. DP (Depth of Coverage): La profondità totale del sequenziamento in quella specifica posizione (quante volte in totale è stata letta quella base).
#     - Note: VAF(Variant Allele Frequency)= NV/DP= percentuale di cellule tumorali che portano quella mutazione
#      
# 3. Stato della mutazione e clonità:
#     1. mutatation_status: Classifica il tipo di alterazione presente. Può essere:
#     - WT (Wild Type): gene sano, non ci sono mutazioni
#     - M (Mutazione puntiforme): mutazione classica, generalmente puntiforme (un singolo nucleotide cambiato, o una piccola inserzione/delezione)
#     - CI_M : NON LO SO - clonal illusion mutations???? mutazioni puntiformi in cui l'algoritmo non è risucito ad assegnare con certezza il timing o l'appartenenza a un clone - incertezza a classificare se sono mutazioni clonali o subclonali
#     - CNA_driver (Copy Number Alteration): gene è diventato un driver del cancro perché è stato fisicamente amplificato (ci sono troppe copie del gene, producendo troppe proteine) o deleto (la cellula ha perso il gene protettivo)
# 
#     2. mult_estimate: La molteplicità della mutazione, ovvero su quante copie fisiche del cromosoma è presente quell'errore genetico (es. 1 o 2)- se la utaizone è avvenuta pre o post duplicaizone
#     3. mutation_call: Indica esplicitamente se la mutazione colpisce un gene "driver" (un gene che guida attivamente lo sviluppo del cancro) riconosciuto dai cataloghi ufficiali (es. PCAWG)
# 
# 4. Timing evolutivo:
#     1. timed: Un flag (Vero/Falso) che indica se è stato possibile datare con successo la comparsa della mutazione
#     2. clock_rank: ordine cronologico, rank del cluster a cui sono assegnati gli eventi tumori avvenuti nello stesso intervallo
#     3. clock_mean / clock_low / clock_high: La stima numerica (media e intervalli di confidenza) del tempo evolutivo in cui la mutazione è emersa.
# 
# 5. Caratteristiche flobali del tumore: 
#     1. class: traiettoria del tumore:
#     - WGD (genoma raddoppiato): il tumore ha preso l'intero genoma (46 cromosomi) e lo ha raddoppiato.
#     - HM (Hypermutated, tumori con un tasso di mutazione altissimo)  
#     - Classic (evoluzione standard).
#     2. is_WGD: Indica se il tumore ha subito un evento di Whole Genome Doubling (raddoppiamento dell'intero genoma).
#     3. ploidy: [3. 2. 4. 5. 6. 1.]: La ploidia totale del tumore, che riflette il numero medio di set di cromosomi (valori alti indicano forte instabilità genomica)- cellula normale è diploide (ploidy=2), cioè ha due copie di ogni cromosoma
#     4. Purity: La purezza del campione, ovvero la percentuale effettiva di cellule tumorali presenti rispetto alle cellule sane infiltrate nel tessuto analizzato.
#     5. best_k:numero di ondate di amplificazioni- il tumore è:
#     - monoclonale: tutte le cellule malate sono geneticamente simili
#     - policlonale: diverse famiglie di cellule tumorali che convivono 
#     6. ttype (Tumor Type, 30): La sigla ufficiale che identifica l'istotipo del cancro (es. LUAD per il polmone, BRCA per il seno, PRAD per la prostata).

# %% [markdown]
# # Analisi esplorativa (EDA)

# %% [markdown]
# ## Missing data e analisi generale

# %%
# Calcolo dei dati mancanti: Numero / Totale
total_rows = len(df)
missing_counts = df.isnull().sum()

missing_stats = pd.DataFrame({
    'Mancanti': missing_counts,
    'Totale_Righe': total_rows,
    'Percentuale (%)': (missing_counts / total_rows) * 100
})

# Mostriamo le colonne ordinate per quelle con più dati mancanti
print(missing_stats.sort_values(by='Mancanti', ascending=False))

# %%
print(df[df['NV'].notnull()]['mutatation_status'].value_counts())
# 1. Conteggio per SEGMENTI (Quante righe totali nel dataset appartengono a ogni categoria)
print("--- CONTEGGIO PER SEGMENTI (RIGHE TOTALI) ---")
status_counts = df['mutatation_status'].value_counts()
print(status_counts)

# Calcoliamo anche la percentuale per avere un'idea più chiara
print("\n--- PERCENTUALE PER SEGMENTI ---")
print((status_counts / len(df) * 100).round(2).astype(str) + ' %')


# 2. Conteggio per PAZIENTI (Quanti pazienti presentano ALMENO UNA volta quel determinato status)
print("\n--- CONTEGGIO PER PAZIENTI UNICI ---")
# Raggruppiamo per status e contiamo i sample_id unici
patient_status_counts = df.groupby('mutatation_status')['sample_id'].nunique().sort_values(ascending=False)
print(patient_status_counts)

# Calcoliamo la percentuale rispetto al numero totale di pazienti
total_patients = df['sample_id'].nunique()
print("\n--- PERCENTUALE SUI PAZIENTI TOTALI ---")
print((patient_status_counts / total_patients * 100).round(2).astype(str) + ' %')

# %% [markdown]
# Le prime colonne identificano le mutazioni puntiformi cioè lo stato M o CI_M: un singolo nucleotide cambiato, o una piccola inserzione/delezione.  
# 
# Lo studio analizza alterazioni strutturali del numero di copie, cioè alterazioni del numero degli alleli nel cromosoma e identifica l'ordine dei segmenti cromosomici. Ogni riga del dataset è un gene che si trova nel segmento alterato.
# 
# I dati fanno riferimento a geni:
# 
# - gene WT (Wild type): gene sano, senza mutazioni
# - alterati da CNA (Copy Number Alteration): il gene è diventato un driver del cancro perché l'intero frammento di cromosoma è stato copiato troppe volte (amplificazione) o cancellato (delezione). Poiché l'intero blocco è anomalo, non c'è una singola "lettera" del DNA mutata di cui puoi calcolare la molteplicità (mult_estimate). L'informazione sul danno strutturale è già contenuta tutta nella colonna karyotype. 
# - gene M o CI_M: un singolo nucleotide (base del DNA= lettera), o una piccola inserzione/delezione 
# Si nota subito che: 
# - i geni sono quasi tutti WT: quando un segmento di cromosoma si duplica o elimina, tutti i geni posizionati in esso vengono duplicati di conseguenza e la maggior parte di essi è sana, sono solo stati trascinati nell'amplificazione del segmento in cui risiedono. 
# 
#   L'algoritmo ha scansionato i geni di ogni segmento alla ricerca delle mutazioni per creare l'ordine cronologico: se sono state duplicate hanno un alta frequenza, se sono nate dopo l'amplificaizone hanno una bassa freq.
#   
#   Tutti i geni all'interno dello stesso segmento ereditano lo stesso clock, anche se sono sani.  
# 
# - Sono presenti molti null in NV, DP, mult_estimate e mutation_call: sono colonne specifiche per mutazioni puntiformi, se la maggior parte dei geni è sana, non c'è una mutazione.
# 
# 

# %%
# ==========================================
# 1. PAZIENTI E TIPI DI TUMORE
# ==========================================
print("-" * 50)
print("1. PAZIENTI E TIPI DI TUMORE (ttype)")
print("-" * 50)

pazienti_totali = df['sample_id'].nunique()
tumori_totali = df['ttype'].nunique()

print(f"Totale Pazienti unici: {pazienti_totali}")
print(f"Totale Tipi di tumore (ttype) diversi: {tumori_totali}")

# Verifichiamo se ogni paziente ha un solo tipo di tumore
tumori_per_paziente = df.groupby('sample_id')['ttype'].nunique()
pazienti_con_piu_tumori = (tumori_per_paziente > 1).sum()

if pazienti_con_piu_tumori == 0:
    print("Conferma: Ogni paziente ha ESATTAMENTE un solo tipo di tumore.")
else:
    print(f"Attenzione: {pazienti_con_piu_tumori} pazienti hanno più di un tipo di tumore.")

# %%
# ==========================================
# 2. GLI EVENTI STRUTTURALI (Isoliamo i Cluster)
# ==========================================
print("\n" + "-" * 50)
print("2. I CLOCK RANK E GLI INTERVALLI DI TEMPO")
print("-" * 50)

# FONDAMENTALE: Rimuoviamo i geni "passeggeri" duplicati. 
# Teniamo solo una riga per ogni cluster (clock_rank) di ogni paziente.
eventi_unici = df.dropna(subset=['clock_rank']).drop_duplicates(subset=['sample_id', 'clock_rank'])

# Quanti clock_rank ci sono in totale?
rank_counts = eventi_unici['clock_rank'].value_counts().sort_index()
print("Frequenza dei Clock Rank (Quante 'ondate' strutturali ci sono):")
for rank, count in rank_counts.items():
    print(f"- Rank {int(rank)}: presente in {count} tumori")

# Che intervalli di tempo coprono mediamente questi rank?
statistiche_tempo = eventi_unici.groupby('clock_rank')['clock_mean'].agg(['min', 'mean', 'max'])
print("\nTempi medi (clock_mean) per ogni Rank (da 0 = nascita a 1 = prelievo):")
print(statistiche_tempo.round(3))

# %% [markdown]
# Si può notare che: 
# - Quasi la totalità della coorte (997 tumori su 1002) presenta almeno una grande "ondata" di amplificazioni strutturali databili (Rank 1).
# 
# - La frequenza crolla drasticamente per i rank successivi- c'è evoluzione saltatoria (saltational evolution / hopeful monsters): la maggior parte dei tumori non accumula danni strutturali in modo continuo e frammentato, ma acquisisce le sue alterazioni cromosomiche principali in uno o pochissimi eventi massicci.

# %%
# ==========================================
# 3. VISUALIZZAZIONE DELLE SOVRAPPOSIZIONI
# ==========================================
# Creiamo un grafico per rispondere alla tua domanda: "quanti intervalli si sovrappongono?"
plt.figure(figsize=(10, 6))
sns.boxplot(data=eventi_unici, x='clock_rank', y='clock_mean', palette='viridis', width=0.6)
sns.stripplot(data=eventi_unici, x='clock_rank', y='clock_mean', color='black', alpha=0.3, size=3)

plt.title('Sovrapposizione dei Tempi Evolutivi (clock_mean) per ogni Rank')
plt.xlabel('Clock Rank (Ordine cronologico)')
plt.ylabel('Tempo Molecolare (clock_mean)')
plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.show()

# %% [markdown]
# Distribuzione del tempo molecolare ($\tau$) per ogni ondata (Cluster-puntini):
# - L'ampiezza del Rank 1: La prima ondata di instabilità (Rank 1) ha una varianza enorme. Può avvenire all'alba della nascita del tumore (vicino a 0.0) o in fasi molto avanzate (vicino a 1.0). 
# - L'evento primario non ha una scadenza fissa.La compressione dei Rank successivi: Come ci si aspetta logicamente, i Rank successivi (2, 3, 4, 5) sono progressivamente schiacciati verso l'alto (verso 1.0). Se un tumore subisce 4 o 5 "terremoti" cromosomici, gli ultimi avvengono per forza in tempi molto recenti, immediatamente prima del prelievo bioptico. Tick Tack ha ordinato cronologicamente gli eventi in modo matematicamente coerente.

# %%
### Check a tumor with k max mean between all patient with that tumor > 1

import pyreadr
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


# Load the rds file
result = pyreadr.read_r('06_Cb_BTM_table.rds')

# pyreadr returns a dictionary; the data is usually under the 'None' key 
# because rds files contain a single object without a specific name.
df = result[None]

df = df[df['mutatation_status']!="WT"]


## remove ttype with few patients
treshold = 10

tumors = df["ttype"].unique()

tumors_and_mean = {}

for t in tumors:
    df_t = df[df["ttype"] == t]
    if len(df_t["sample_id"].unique())<treshold:
        continue
    df_t = df_t[["sample_id","clock_rank"]]
    df_t_p = df_t.groupby("sample_id")["clock_rank"].max()
    tumors_and_mean[t] = df_t_p.mean()
    
max_key = max(tumors_and_mean, key=tumors_and_mean.get)

print(max_key, tumors_and_mean[max_key])

    


