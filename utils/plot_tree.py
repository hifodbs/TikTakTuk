import graphviz
import os
def print_tree(t, name):
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
    
    
    
def print_tree_second(t, name):
    # Remove the "produced/" prefix from here, so it will be controlled it from the main script
    if name.endswith(".gv"):
        filename = name[:-3]
    else:
        filename = name

    d = graphviz.Digraph(filename=filename, format="pdf")

    # layout
    d.attr(rankdir="TB")         # from L to R tree
    d.attr(ranksep="2.5")        # Distance bt ranks
    d.attr(nodesep="0.4")        #distance bt nodes
       
    # color and shape of the cell
    d.attr("node", shape="box", style="rounded,filled", fillcolor="aliceblue", 
           fontname="Helvetica", fontsize="11", margin="0.2,0.1")

    for k, v in t.items():
        # Graphviz aligns the cell without the box
        with d.subgraph(name=f"rank_{k}") as s:
            s.attr(rank="same")
            for label_name in v:
                s.node(label_name, label=label_name)

    # Add edges
    d.attr("edge", color="gray75", penwidth="0.8", arrowsize="0.6")
    
    max_rank = max(t.keys())
    for r in range(1, max_rank):
        if r in t and (r + 1) in t:
            for n1 in t[r]:
                for n2 in t[r + 1]:
                    d.edge(n1, n2)

    # Save the file (which automatically generates the .gv and .pdf)
    d.render(cleanup=False) 
    d.view
    print(f"Graph rendered: {filename}.gv and {filename}.pdf")