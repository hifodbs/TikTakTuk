from pyvis.network import Network
import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt




def plot(G):
    # 1. Initialize PyVis network (directed)
    net = Network(height="750px", width="100%", directed=True, notebook=True)

    # 2. Convert your NetworkX graph G to PyVis
    net.from_nx(G)

    # 3. Force Hierarchical layout & enable physics configuration UI
    net.set_options("""
    var options = {
    "layout": {
        "hierarchical": {
        "enabled": true,
        "direction": "UD",
        "sortMethod": "directed"
        }
    },
    "physics": {
        "hierarchicalRepulsion": {
        "nodeDistance": 150
        },
        "solver": "hierarchicalRepulsion"
    }
    }
    """)

    # 4. Generate and save the interactive HTML file
    net.show("dag_visualization.html")
    
def plot_dot(G):
    # You must have Graphviz installed on your system (e.g., `brew install graphviz` or `sudo apt install graphviz`)
    # Convert NetworkX to a Graphviz AGraph
    A = nx.drawing.nx_agraph.to_agraph(G)

    # Apply hierarchical layout configuration
    A.layout(prog='dot') 

    # Save as SVG (scalable and crisp) or PDF
    A.draw('my_hierarchy.svg')
    
def plot_old_style(gad):
    for layer, nodes in enumerate(nx.topological_generations(gad)):
    #     # `multipartite_layout` expects the layer as a node attribute, so add the
    #     # numeric layer value as a node attribute
        for node in nodes:
            gad.nodes[node]["layer"] = layer

    # Compute the multipartite_layout using the "layer" node attribute
    pos = nx.multipartite_layout(gad, subset_key="layer")

    fig, ax = plt.subplots(figsize=(12, 8))
    nx.draw_networkx(gad, pos=pos, ax=ax, node_color='lightblue', node_size=1500, font_size=10)
    ax.set_title("DAG layout in topological order")
    fig.tight_layout()
    
    # # --- SALVATAGGIO NELLA CARTELLA SPECIFICA ---
    # # 1. Definisci il nome della cartella
    # output_dir = "produced_all_process"
    
    # # 2. Crea la cartella se non esiste già
    # if not os.path.exists(output_dir):
    #     os.makedirs(output_dir)
        
    # # 3. Definisci il nome del file e unisci il percorso della cartella
    # nome_tumore = tumore if tumore is not None else "All"
    # file_immagine = os.path.join(output_dir, f"DAG_evoluzione_{nome_tumore}.png")
    
    # # 4. Salva fisicamente il file
    # plt.savefig(file_immagine, dpi=300, bbox_inches='tight')
    # print(f"--> Grafo salvato con successo nel percorso: {file_immagine}")
    
    # # 5. Mostra l'immagine a schermo (opzionale)
    plt.show()
    
def export_data(gad, tumore):
    df_name = pd.DataFrame({"gene":list(nx.topological_sort(gad))})
    df_name.to_csv(f"file_for_fluorish/genes {tumore}.csv", header=False, index=False)
    
    df_edges = pd.DataFrame(list(gad.edges(data="weight")), columns=["source", "target", "weight"])
    df_edges.to_csv(f"file_for_fluorish/edges {tumore}.csv", header=False, index=False)
    
    

import dash
from dash import html, Input, Output
import dash_cytoscape as cyto
import networkx as nx

# Carica i layout aggiuntivi (incluso dagre)
cyto.load_extra_layouts()

def create_interactive_dag_app(G, node_label_attr=None, edge_weight_attr='weight'):
    """
    Crea un'applicazione Dash interattiva per un DAG.
    - Layout compatto 'dagre'
    - Al click su un nodo evidenzia solo i successori (nodi a valle) e resetta cliccando lo sfondo.
    """
    
    # 1. Prepariamo gli elementi per Cytoscape
    cyto_elements = []
    
    # Nodi
    for node, data in G.nodes(data=True):
        label = str(data.get(node_label_attr)) if node_label_attr else str(node)
        cyto_elements.append({
            'data': {
                'id': str(node), 
                'label': label
            }
        })

    # Archi
    for source, target, data in G.edges(data=True):
        weight = data.get(edge_weight_attr, 1.0)
        cyto_elements.append({
            'data': {
                'id': f"{source}-{target}",
                'source': str(source), 
                'target': str(target),
                'weight': float(weight)
            }
        })

    # 2. Inizializziamo l'applicazione Dash
    app = dash.Dash(__name__)

    # Foglio di stile base
    default_stylesheet = [
        {
            'selector': 'node',
            'style': {
                'label': 'data(label)',
                'background-color': '#0074D9',
                'color': '#2D3748',
                'font-size': '250px',
                'font-weight': 'bold',
                'font-family': 'sans-serif',
                'width': '250px',
                'height': '250px',
                'text-halign': 'left',      # Aligns the text horizontally to the left of the node
                'text-valign': 'center',    # Aligns the text vertically with the center of the node
                'text-margin-x': -8,        # Pushes the text slightly further left so it doesn't overlap the node border
                'opacity': 1.0,
                'transition-property': 'opacity, background-color',
                'transition-duration': '0.15s'
            }
        },
        {
            'selector': 'edge',
            'style': {
                'width': 'mapData(weight, 0, 10, 1.5, 6)', 
                'line-color': '#A0AEC0',
                'target-arrow-color': '#A0AEC0',
                'target-arrow-shape': 'triangle',
                'curve-style': 'bezier',
                'opacity': 0.6,
                'transition-property': 'opacity, line-color, target-arrow-color',
                'transition-duration': '0.15s'
            }
        }
    ]

    # Layout dell'app
    app.layout = html.Div(style={'fontFamily': 'sans-serif', 'padding': '20px'}, children=[
        html.H2("Progressione Clonale - Analisi delle Traiettorie A Valle", style={'textAlign': 'center', 'color': '#2D3748'}),
        html.P("Clicca su un gene per evidenziare solo i suoi successori (discendenti). Clicca sullo sfondo (o seleziona un altro nodo) per resettare.", 
               style={'textAlign': 'center', 'color': '#718096', 'fontSize': '14px'}),
        
        html.Div(style={'display': 'flex', 'justifyContent': 'center'}, children=[
            cyto.Cytoscape(
                id='cytoscape-dag',
                elements=cyto_elements,
                style={'width': '95%', 'height': '800px', 'border': '1px solid #E2E8F0', 'borderRadius': '12px', 'backgroundColor': '#F7FAFC'},
                
                # MODIFICA QUI IL LAYOUT:
                layout={
                    'name': 'dagre',
                    'nodeDimensionsIncludeLabels': False,
                    'animate': True,
                    # rankSep imposta la distanza verticale tra i layer (es. 150-200 pixel)
                    'rankSep': 2000,  
                    # nodeSep mantiene i nodi vicini tra loro sullo stesso livello (es. 30-50 pixel)
                    'nodeSep': 50,   
                    # rankDir definisce la direzione del flusso: 'TB' (dall'alto in basso) o 'LR' (da sinistra a destra)
                    'rankDir': 'LR',
                    'ranker': 'network-simplex',
                },
                stylesheet=default_stylesheet,
                autounselectify=False,
                boxSelectionEnabled=False
            )
        ])
    ])

    # 3. Callback per gestire l'evidenziazione e il reset istantaneo dallo sfondo
    @app.callback(
        Output('cytoscape-dag', 'stylesheet'),
        Input('cytoscape-dag', 'selectedNodeData') # <-- Cambiato qui!
    )
    def generate_stylesheet(selected_nodes):
        # Se non c'è nessun nodo selezionato (es. cliccando sullo sfondo o deselezionando), resettiamo
        if not selected_nodes:
            return default_stylesheet

        # Prendiamo l'ultimo nodo effettivamente cliccato/selezionato
        clicked_node_id = selected_nodes[-1]['id']
        
        # Troviamo i successori e i predecessori nel grafo NetworkX
        connected_nodes = set([clicked_node_id])
        if G.has_node(clicked_node_id):
            connected_nodes.update(G.successors(clicked_node_id))
            connected_nodes.update(G.predecessors(clicked_node_id))

        # Creiamo un foglio di stile dinamico basato sul click
        dynamic_rules = [
            # 1. Di default, "spegniamo" tutti i nodi e gli archi
            {
                'selector': 'node',
                'style': {'opacity': 0.15}
            },
            {
                'selector': 'edge',
                'style': {'opacity': 0.05}
            },
            # 2. Riaccendiamo e coloriamo il nodo cliccato e quelli connessi
            {
                'selector': ','.join([f'node[id = "{n}"]' for n in connected_nodes]),
                'style': {
                    'opacity': 1.0,
                    'background-color': '#E53E3E', # Rosso vivo per i nodi attivi
                    'color': '#2D3748'
                }
            },
            # Evidenziamo nello specifico il nodo cliccato rispetto ai vicini
            {
                'selector': f'node[id = "{clicked_node_id}"]',
                'style': {
                    'background-color': '#319795', # Petrolio/Teal per il "centro" del click
                    'width': '270px',              # Adattato ai nodi giganti da 250px del tuo foglio di stile
                    'height': '270px'
                }
            },
            # 3. Riaccendiamo solo gli archi che collegano direttamente il nodo cliccato ai suoi vicini
            {
                'selector': ','.join([
                    f'edge[source = "{clicked_node_id}"], edge[target = "{clicked_node_id}"]'
                ]),
                'style': {
                    'opacity': 1.0,
                    'line-color': '#E53E3E',
                    'target-arrow-color': '#E53E3E',
                    'width': 3
                }
            }
        ]

        # Uniamo lo stile base con le regole dinamiche
        return default_stylesheet + dynamic_rules
    return app