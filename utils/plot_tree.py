import graphviz
import os
def print_tree(t, name):
    
    os.makedirs("./produced", exist_ok=True) 
    # Initialize graph
    d = graphviz.Digraph(filename="./produced/"+name)

    # 1. FIX SPACING: Spread layers vertically and nodes horizontally
    d.attr(ranksep="2.5")  # Distance between layers
    d.attr(nodesep="0.2")  # Distance between nodes on the same layer
    d.attr(rankdir="TB")  # Strict Top-to-Bottom layout
    

    # 2. FIX TEXT SIZE: Force a clean box shape and stable font size
    d.attr("node", shape="box", fontsize="14", width="0.3", height="0.5")

    layers = t

    # Create nodes within their respective layers
    for k, v in layers.items():
        with d.subgraph(name=f"cluster_{k}") as s:
            s.attr(rank="same")
            s.attr(label=f"Rank {k}", fontsize="16", color="lightgrey")
            for gene in v:
                if len(v)>20:
                    # Splits "TP53" into "T\nP\n5\n3"
                    label_name = "\n".join(list(gene))
                else:
                    label_name = gene
                
                # Pass vertical_name as the label, but keep 'gene' as the internal ID
                s.node(gene, label=label_name)

    # Add edges
    max_rank = max(layers)

    for r in range(1, max_rank):
        for n1 in layers[r]:
            for n2 in layers[r + 1]:
                # NOTE: If your data has explicit parent-child mappings,
                # replace this nested loop with those explicit pairs to avoid the hairball effect.
                d.edge(n1, n2)


    d = d.unflatten(stagger=3)  
    d.view()
    



#helper fcs

def generate_simple_sequential_edges(clusters_dict):
    """Genera archi 'tutti-a-tutti' tra il Rank 1 e il Rank 2, il Rank 2 e il Rank 3, ecc."""
    edges = []
    ranks = sorted(list(clusters_dict.keys()))
    
    for i in range(len(ranks) - 1):
        r1 = ranks[i]
        r2 = ranks[i + 1]
        
        for c1 in clusters_dict[r1].keys():
            for c2 in clusters_dict[r2].keys():
                edges.append({
                    'source': (r1, c1),
                    'target': (r2, c2)
                })
    return edges


def generate_rank_edges(input_dict):
    """Genera archi univoci che collegano l'intero Rank 1 all'intero Rank 2, ecc."""
    edges = []
    ranks = sorted(list(input_dict.keys()))
    
    for i in range(len(ranks) - 1):
        edges.append({
            'source_rank': ranks[i],
            'target_rank': ranks[i + 1]
        })
    return edges

# ==========================================
# FUNZIONE DI CENTRAMENTO ASSOLUTO: grafico
# ==========================================
def _get_middle_node_visual(clusters_dict, rank, is_gene=False):
    """
    Trova il nodo fisicamente centrale leggendo i dati 
    nell'esatto ordine in cui verranno disegnati a schermo, 
    senza ordinamenti alfabetici che sfaserebbero il grafico.
    """
    nodes = []
    # Usiamo l'ordine naturale (insertion order) del dizionario
    for c_id, elements in clusters_dict.get(rank, {}).items():
        if is_gene:
            # Per i geni, elements è un dict: {gene: [pazienti]}
            for gene in elements.keys():
                nodes.append(f"{gene}_R{rank}_C{c_id}")
        else:
            # Per i pazienti, elements è una lista di numeri
            for p in elements:
                nodes.append(f"P_{p}_R{rank}_C{c_id}")
    
    if nodes:
        # Restituisce l'elemento esattamente a metà della riga visiva
        return nodes[len(nodes) // 2]
    return None

# ==========================================
# PLOT PAZIENTI
# ==========================================
def plot_patients_clusters_simple(patient_clusters, edges, name):
    os.makedirs("./produced_all_process", exist_ok=True) 
    d = graphviz.Digraph(filename="./produced_all_process/"+name)

    d.attr(compound="true") 
    d.attr(ranksep="2.5", nodesep="0.3", rankdir="TB")  
    d.attr("node", shape="box", style="filled", fillcolor="lightyellow", fontsize="14")

    # 1. Disegna Nodi e Cluster
    for rank, clusters in patient_clusters.items():
        with d.subgraph(name=f"cluster_rank_{rank}") as s_rank:
            s_rank.attr(rank="same")
            s_rank.attr(label=f"Rank {rank}", fontsize="20", style="dashed", color="grey", labeljust="c")

            for c_id, patients in clusters.items():
                if not patients:
                    continue
                with s_rank.subgraph(name=f"cluster_R{rank}_C{c_id}") as s_cluster:
                    s_cluster.attr(label=f"Cluster {c_id}", style="filled", color="lightgrey")
                    for p_num in patients:
                        node_id = f"P_{p_num}_R{rank}_C{c_id}"
                        s_cluster.node(node_id, label=f"Paziente {p_num}")

    # 2. Crea la freccia perfettamente centrale
    for edge in edges:
        r1 = edge['source_rank']
        r2 = edge['target_rank']
        
        n1 = _get_middle_node_visual(patient_clusters, r1, is_gene=False)
        n2 = _get_middle_node_visual(patient_clusters, r2, is_gene=False)
        
        if n1 and n2:
            # weight="100" è il segreto che trasforma la freccia in una colonna portante verticale
            d.edge(n1, n2, 
                   ltail=f"cluster_rank_{r1}", 
                   lhead=f"cluster_rank_{r2}", 
                   color="#222222", 
                   penwidth="2.5", 
                   weight="100")

    d.view()

# ==========================================
# PLOT GENI
# ==========================================
def plot_genes_clusters_simple(gene_clusters, edges, name):
    os.makedirs("./produced_all_process", exist_ok=True) 
    d = graphviz.Digraph(filename="./produced_all_process/"+name)

    d.attr(compound="true") 
    d.attr(ranksep="2.5", nodesep="0.3", rankdir="TB")  
    d.attr("node", shape="box", style="filled", fillcolor="lightblue", fontsize="14")

    # 1. Disegna Nodi e Cluster
    for rank, clusters in gene_clusters.items():
        with d.subgraph(name=f"cluster_rank_{rank}") as s_rank:
            s_rank.attr(rank="same")
            s_rank.attr(label=f"Rank {rank}", fontsize="20", style="dashed", color="grey", labeljust="c")

            for c_id, genes_dict in clusters.items():
                if not genes_dict:
                    continue
                with s_rank.subgraph(name=f"cluster_R{rank}_C{c_id}") as s_cluster:
                    s_cluster.attr(label=f"Cluster {c_id}", style="filled", color="lightgrey")
                    for gene, p_list in genes_dict.items():
                        node_id = f"{gene}_R{rank}_C{c_id}"
                        p_string = ", ".join(map(str, sorted(p_list, key=int)))
                        label_text = f"{gene}\n({p_string})"
                        s_cluster.node(node_id, label=label_text)

    # 2. Crea la freccia perfettamente centrale
    for edge in edges:
        r1 = edge['source_rank']
        r2 = edge['target_rank']
        
        n1 = _get_middle_node_visual(gene_clusters, r1, is_gene=True)
        n2 = _get_middle_node_visual(gene_clusters, r2, is_gene=True)
        
        if n1 and n2:
            d.edge(n1, n2, 
                   ltail=f"cluster_rank_{r1}", 
                   lhead=f"cluster_rank_{r2}", 
                   color="#222222", 
                   penwidth="2.5", 
                   weight="100")

    d.view()
    
    
# PMI
    
def _get_first_gene_node_in_cluster(genes_dict, rank, c_id):
    """Funzione di supporto per trovare il primo nodo utile in un cluster specifico."""
    if genes_dict:
        first_gene = list(genes_dict.keys())[0]
        return f"{first_gene}_R{rank}_C{c_id}"
    return None

def plot_genes_clusters_pmi(gene_clusters, pmi_edges, name):
    import os
    import graphviz
    
    os.makedirs("./produced_all_process", exist_ok=True) 
    d = graphviz.Digraph(filename="./produced_all_process/"+name)

    d.attr(compound="true") 
    d.attr(ranksep="2.5", nodesep="0.3", rankdir="TB")  
    d.attr("node", shape="box", style="filled", fillcolor="lightblue", fontsize="14")

    # Dizionario per tenere traccia delle ancore di ogni Rank
    rank_anchors = {}

    # 1. Disegna tutti i Cluster con i geni all'interno
    for rank, clusters in gene_clusters.items():
        with d.subgraph(name=f"cluster_rank_{rank}") as s_rank:
            s_rank.attr(rank="same")
            s_rank.attr(label=f"Rank {rank}", fontsize="20", style="dashed", color="grey", labeljust="c")

            # CREAZIONE ANCORA: Un punto minuscolo e completamente invisibile
            anchor_id = f"anchor_R{rank}"
            s_rank.node(anchor_id, style="invis", shape="point", width="0", height="0")
            rank_anchors[rank] = anchor_id

            for c_id, genes_dict in clusters.items():
                if not genes_dict:
                    continue
                with s_rank.subgraph(name=f"cluster_R{rank}_C{c_id}") as s_cluster:
                    s_cluster.attr(label=f"Cluster {c_id}", style="filled", color="lightgrey")
                    for gene, p_list in genes_dict.items():
                        node_id = f"{gene}_R{rank}_C{c_id}"
                        p_string = ", ".join(map(str, sorted(p_list, key=int)))
                        label_text = f"{gene}\n({p_string})"
                        s_cluster.node(node_id, label=label_text)

    # COLLONNA VERTEBRALE INVISIBILE: Forza l'allineamento verticale dei Rank
    sorted_ranks = sorted(rank_anchors.keys())
    for i in range(len(sorted_ranks) - 1):
        d.edge(rank_anchors[sorted_ranks[i]], rank_anchors[sorted_ranks[i+1]], style="invis")

    # 2. Crea le frecce evolutive tra cluster basate sulla PMI
    for edge in pmi_edges:
        r_s, c_s = edge['source']
        r_t, c_t = edge['target']
        w = edge['weight']
        pmi = edge['pmi']
        
        n1 = _get_first_gene_node_in_cluster(gene_clusters.get(r_s, {}).get(c_s), r_s, c_s)
        n2 = _get_first_gene_node_in_cluster(gene_clusters.get(r_t, {}).get(c_t), r_t, c_t)
        
        if n1 and n2:
            label_text = f"{w} pts\n PMI: {pmi:.1f}"  # 
            thickness = str(1.0 + (pmi * 1.5))
            
            d.edge(n1, n2, 
                   ltail=f"cluster_R{r_s}_C{c_s}", 
                   lhead=f"cluster_R{r_t}_C{c_t}", 
                   label=label_text,
                   fontsize="11",
                   color="#333333", 
                   penwidth=thickness)

    d.view()
    
    
#GLOBAL 
    
def plot_genes_clusters_pmi_global(gene_clusters, pmi_edges, name):
    os.makedirs("./nuova_cartella", exist_ok=True) 
    d = graphviz.Digraph(filename="./nuova_cartella/"+name)

    d.attr(compound="true") 
    d.attr(ranksep="2.5", nodesep="0.3", rankdir="TB")  
    d.attr("node", shape="box", style="filled", fillcolor="lightblue", fontsize="14")

    # 1. Disegna i Cluster
    for rank, clusters in gene_clusters.items():
        with d.subgraph(name=f"cluster_rank_{rank}") as s_rank:
            s_rank.attr(rank="same")
            s_rank.attr(label=f"Rank {rank}", fontsize="20", style="dashed", color="grey", labeljust="c")
            for c_id, genes_dict in clusters.items():
                if not genes_dict: continue
                with s_rank.subgraph(name=f"cluster_R{rank}_C{c_id}") as s_cluster:
                    s_cluster.attr(label=f"Cluster {c_id}", style="filled", color="lightgrey")
                    for gene, p_list in genes_dict.items():
                        node_id = f"{gene}_R{rank}_C{c_id}"
                        # Etichetta pulita: solo il gene e il numero di pazienti che lo hanno
                        p_string = ", ".join(map(str, sorted(p_list, key=int)))
                        s_cluster.node(node_id, label=f"{gene}\n({p_string})")

    # 2. Crea le frecce basate sulla PMI globale
    for edge in pmi_edges:
        r_s, c_s = edge['source']
        r_t, c_t = edge['target']
        pmi = edge['pmi'] # PMI è la forza della correlazione
        
        n1 = _get_first_gene_node_in_cluster(gene_clusters.get(r_s, {}).get(c_s), r_s, c_s)
        n2 = _get_first_gene_node_in_cluster(gene_clusters.get(r_t, {}).get(c_t), r_t, c_t)
        
        if n1 and n2:
            # ETICHETTA: Mostriamo solo il valore della correlazione (PMI)
            label_text = f"PMI: {pmi:.2f}"
            # SPESSORE: Più alta è la PMI, più spessa è la freccia
            thickness = str(max(1.0, min(5.0, pmi * 2)))
            
            d.edge(n1, n2, 
                   ltail=f"cluster_R{r_s}_C{c_s}", 
                   lhead=f"cluster_R{r_t}_C{c_t}", 
                   label=label_text,
                   fontsize="10",
                   color="#333333", 
                   penwidth=thickness)

    d.view()