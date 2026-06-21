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
    
    
import graphviz

