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
    A.layout(prog="dot")

    # Save as SVG (scalable and crisp) or PDF
    A.draw("my_hierarchy.svg")


def plot_old_style(gad):
    for layer, nodes in enumerate(nx.topological_generations(gad)):
        #     # `multipartite_layout` expects the layer as a node attribute, so add the
        #     # numeric layer value as a node attribute
        for node in nodes:
            gad.nodes[node]["layer"] = layer

    # Compute the multipartite_layout using the "layer" node attribute
    pos = nx.multipartite_layout(gad, subset_key="layer")

    fig, ax = plt.subplots(figsize=(12, 8))
    nx.draw_networkx(
        gad, pos=pos, ax=ax, node_color="lightblue", node_size=1500, font_size=10
    )
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
    df_name = pd.DataFrame({"gene": list(nx.topological_sort(gad))})
    df_name.to_csv(f"file_for_fluorish/genes {tumore}.csv", header=False, index=False)

    df_edges = pd.DataFrame(
        list(gad.edges(data="weight")), columns=["source", "target", "weight"]
    )
    df_edges.to_csv(f"file_for_fluorish/edges {tumore}.csv", header=False, index=False)


import dash
from dash import html, Input, Output
import dash_cytoscape as cyto
import networkx as nx

# Carica i layout aggiuntivi (incluso dagre)
cyto.load_extra_layouts()


def create_interactive_dag_app(G, node_label_attr=None, edge_weight_attr="weight"):
    """
    Crea un'applicazione Dash interattiva per un DAG.
    - Layout compatto 'dagre'
    - Al click su un nodo evidenzia solo i successori (nodi a valle) e resetta cliccando lo sfondo.
    """

    node_size = "250px"
    text_size = "250px"
    space_between_nodes = 50
    space_between_columns = 2000
    color_node = "#0074D9"
    color_parent = "#9735C1"
    color_child = "#D92929"
    color_node_selected = "#319795"
    size_edge_hilight = 10

    # 1. Prepariamo gli elementi per Cytoscape
    cyto_elements = []

    # Nodi
    for node, data in G.nodes(data=True):
        label = str(data.get(node_label_attr)) if node_label_attr else str(node)
        cyto_elements.append({"data": {"id": str(node), "label": label}})

    # Archi
    for source, target, data in G.edges(data=True):
        weight = data.get(edge_weight_attr, 1.0)
        cyto_elements.append(
            {
                "data": {
                    "id": f"{source}-{target}",
                    "source": str(source),
                    "target": str(target),
                    "weight": float(weight),
                }
            }
        )

    # 2. Inizializziamo l'applicazione Dash
    app = dash.Dash(__name__)

    # Foglio di stile base
    default_stylesheet = [
        {
            "selector": "node",
            "style": {
                "label": "data(label)",
                "background-color": color_node,
                "color": "#2D3748",
                "font-size": text_size,
                "font-weight": "bold",
                "font-family": "sans-serif",
                "width": node_size,
                "height": node_size,
                "text-halign": "left",
                "text-valign": "center",
                "text-margin-x": -8,
                # --- SFONDO DEL TESTO (LEGENDA/ETICHETTA) ---
                "text-background-color": "#FFFFFF",  # Colore di sfondo bianco
                "text-background-opacity": 0.75,  # Opacità (0.75 = 75% opaco, lascia intravedere il sottofondo)
                "text-background-shape": "roundrectangle",  # Bordi leggermente arrotondati per un look più pulito
                "text-background-padding": "3px",  # Crea un piccolo cuscinetto di spazio attorno alle lettere
                # --------------------------------------------
                "opacity": 1.0,
                "transition-property": "opacity, background-color",
                "transition-duration": "0.15s",
            },
        },
        {
            "selector": "edge",
            "style": {
                "width": "mapData(weight, 0, 10, 1.5, 6)",
                "line-color": "#A0AEC0",
                "target-arrow-color": "#A0AEC0",
                "target-arrow-shape": "triangle",
                "curve-style": "bezier",
                "opacity": 0.6,
                "transition-property": "opacity, line-color, target-arrow-color",
                "transition-duration": "0.15s",
            },
        },
    ]

    # Layout dell'app
    app.layout = html.Div(
        style={"fontFamily": "sans-serif", "padding": "20px"},
        children=[
            html.H2(
                "Clonal progression — trajectory analysis",
                style={"textAlign": "center", "color": "#2D3748"},
            ),
            html.P(
                "Click on a gene to highlight its direct parents and direct children.",
                style={"textAlign": "center", "color": "#718096", "fontSize": "14px"},
            ),
            html.Div(
                style={
                    "display": "flex",
                    "justifyContent": "center",
                    "position": "relative",
                },
                children=[
                    # --- LEGENDA DENTRO IL GRAFO (TOP LEFT) ---
                    html.Div(
                        style={
                            "position": "absolute",
                            "top": "20px",
                            "left": "60px",
                            "zIndex": "999",  # Mantiene la legenda sopra il canvas di Cytoscape
                            "backgroundColor": "rgba(255, 255, 255, 0.9)",  # Sfondo leggermente trasparente
                            "padding": "12px 16px",
                            "borderRadius": "8px",
                            "border": "1px solid #E2E8F0",
                            "boxShadow": "0 4px 6px -1px rgba(0, 0, 0, 0.1)",
                            "display": "flex",
                            "flexDirection": "column",
                            "gap": "8px",
                            "fontSize": "12px",
                            "fontWeight": "bold",
                            "fontFamily": "sans-serif",
                        },
                        children=[
                            html.Div(
                                style={"display": "flex", "alignItems": "center"},
                                children=[
                                    html.Span(
                                        style={
                                            "display": "inline-block",
                                            "width": "10px",
                                            "height": "10px",
                                            "backgroundColor": color_node_selected,
                                            "borderRadius": "50%",
                                            "marginRight": "8px",
                                        }
                                    ),
                                    html.Span(
                                        "Selected node", style={"color": "#2D3748"}
                                    ),
                                ],
                            ),
                            html.Div(
                                style={"display": "flex", "alignItems": "center"},
                                children=[
                                    html.Span(
                                        style={
                                            "display": "inline-block",
                                            "width": "10px",
                                            "height": "10px",
                                            "backgroundColor": color_parent,
                                            "borderRadius": "50%",
                                            "marginRight": "8px",
                                        }
                                    ),
                                    html.Span(
                                        "Direct parents (left)",
                                        style={"color": "#2D3748"},
                                    ),
                                ],
                            ),
                            html.Div(
                                style={"display": "flex", "alignItems": "center"},
                                children=[
                                    html.Span(
                                        style={
                                            "display": "inline-block",
                                            "width": "10px",
                                            "height": "10px",
                                            "backgroundColor": color_child,
                                            "borderRadius": "50%",
                                            "marginRight": "8px",
                                        }
                                    ),
                                    html.Span(
                                        "Direct children (right)",
                                        style={"color": "#2D3748"},
                                    ),
                                ],
                            ),
                            html.Div(
                                style={"display": "flex", "alignItems": "center"},
                                children=[
                                    html.Span(
                                        style={
                                            "display": "inline-block",
                                            "width": "10px",
                                            "height": "10px",
                                            "backgroundColor": color_node,
                                            "borderRadius": "50%",
                                            "marginRight": "8px",
                                        }
                                    ),
                                    html.Span("Not related nodes", style={"color": "#2D3748"}),
                                ],
                            ),
                        ],
                    ),
                    cyto.Cytoscape(
                        id="cytoscape-dag",
                        elements=cyto_elements,
                        style={
                            "width": "95%",
                            "height": "800px",
                            "border": "1px solid #E2E8F0",
                            "borderRadius": "12px",
                            "backgroundColor": "#F7FAFC",
                        },
                        # MODIFICA QUI IL LAYOUT:
                        layout={
                            "name": "dagre",
                            "nodeDimensionsIncludeLabels": False,
                            "animate": True,
                            # rankSep imposta la distanza verticale tra i layer (es. 150-200 pixel)
                            "rankSep": space_between_columns,
                            # nodeSep mantiene i nodi vicini tra loro sullo stesso livello (es. 30-50 pixel)
                            "nodeSep": space_between_nodes,
                            # rankDir definisce la direzione del flusso: 'TB' (dall'alto in basso) o 'LR' (da sinistra a destra)
                            "rankDir": "LR",
                            "ranker": "network-simplex",
                        },
                        stylesheet=default_stylesheet,
                        autounselectify=False,
                        boxSelectionEnabled=False,
                    ),
                ],
            ),
        ],
    )

    # 3. Callback per gestire l'evidenziazione e il reset istantaneo dallo sfondo
    @app.callback(
        Output("cytoscape-dag", "stylesheet"),
        Input("cytoscape-dag", "selectedNodeData"),
    )
    def update_styles(selected_nodes):
        # Se non c'è nessun nodo selezionato, mostriamo lo stile di default
        if not selected_nodes:
            return default_stylesheet

        clicked_node_id = selected_nodes[-1]["id"]

        # 1. Trova i genitori diretti (immediate predecessors)
        direct_parents = set()
        if G.has_node(clicked_node_id):
            direct_parents.update(G.predecessors(clicked_node_id))

        # 2. Trova solo i figli diretti (immediate successors)
        direct_children = set()
        if G.has_node(clicked_node_id):
            direct_children.update(G.successors(clicked_node_id))

        # Inizializziamo le regole dinamiche spegnendo il resto del grafo di default
        dynamic_rules = [
            {"selector": "node", "style": {"opacity": 0.15}},
            {"selector": "edge", "style": {"opacity": 0.05}},
            # Il nodo cliccato (sorgente) diventa verde petrolio (sempre presente)
            {
                "selector": f'node[id = "{clicked_node_id}"]',
                "style": {
                    "opacity": 1.0,
                    "background-color": color_node_selected,
                    "width": node_size,
                    "height": node_size,
                },
            },
        ]

        # --- COLORAZIONE NODI E ARCHI CONDIZIONALE ---

        # A. Aggiungi figli diretti solo se esistono (nodi non-foglia)
        if direct_children:
            dynamic_rules.append(
                {
                    "selector": ",".join(
                        [f'node[id = "{n}"]' for n in direct_children]
                    ),
                    "style": {
                        "opacity": 1.0,
                        "background-color": color_child,
                        "color": "#2D3748",
                    },
                }
            )

            child_edges = [
                f'edge[id = "{clicked_node_id}-{c}"]' for c in direct_children
            ]
            dynamic_rules.append(
                {
                    "selector": ",".join(child_edges),
                    "style": {
                        "opacity": 1.0,
                        "line-color": color_child,
                        "target-arrow-color": color_child,
                        "width": size_edge_hilight,
                    },
                }
            )

        # B. Aggiungi genitori diretti solo se esistono (nodi non-radice)
        if direct_parents:
            dynamic_rules.append(
                {
                    "selector": ",".join([f'node[id = "{n}"]' for n in direct_parents]),
                    "style": {
                        "opacity": 1.0,
                        "background-color": color_parent,
                        "color": "#2D3748",
                    },
                }
            )

            parent_edges = [
                f'edge[id = "{p}-{clicked_node_id}"]' for p in direct_parents
            ]
            dynamic_rules.append(
                {
                    "selector": ",".join(parent_edges),
                    "style": {
                        "opacity": 1.0,
                        "line-color": color_parent,
                        "target-arrow-color": color_parent,
                        "width": size_edge_hilight,
                    },
                }
            )

        return default_stylesheet + dynamic_rules

    return app
