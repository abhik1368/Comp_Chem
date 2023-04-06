import networkx as nx
from networkx.readwrite import cytoscape_data
from rdkit import Chem
from rdkit.Chem import RDConfig
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import Descriptors
from scipy.spatial.distance import euclidean
import numpy as np
from urllib import parse
import ipycytoscape

def smi2svg(smi):
    mol = Chem.MolFromSmiles(smi)
    try:
        Chem.rdmolops.Kekulize(mol)
    except:
        pass
    drawer = rdMolDraw2D.MolDraw2DSVG(500, 300)
    AllChem.Compute2DCoords(mol)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText().replace("svg:", "")
    return svg

def smi2image(smi):
    svg_string = smi2svg(smi)
    impath = 'data:image/svg+xml;charset=utf-8,' + parse.quote(svg_string, safe="")
    return impath

# load SMILES file 
# format: 1st column SMILES | 2nd column compound name
smis, nams, label = [], [] ,[]
with open("sample.smi.txt", 'r') as smi_nam_list:
    for smi_nam in smi_nam_list:
        smis.append(smi_nam.rstrip().split('\t')[0])
        nams.append(smi_nam.rstrip().split('\t')[1]) 
        label.append(smi_nam.rstrip().split('\t')[2]) 
print(label)
print(smis)
print(str(len(smis))+" molecules loaded from SMILES file")
g = nx.Graph()

# add node
for smi, nam, label in zip(smis, nams, label):
    g.add_node(smi, img=smi2image(smi), name=nam, label=label)
    


# Get physchem descriptors and calculate Euclidean distance
mols = [Chem.MolFromSmiles(x) for x in smis]
physchem_descs = [[desc(mol) for _, desc in Descriptors.descList] for mol in mols]



# Calculate all pairwise distances
distances = []
for i in range(len(smis)):
    for j in range(i):
        dist = euclidean(physchem_descs[i], physchem_descs[j])
        distances.append(dist)

# Set the threshold as the 75th percentile of distances
threshold = np.percentile(distances, 90)

d_min, d_max = min(distances), max(distances)

# Add edges with distances as edge weights and labels
for i in range(len(smis)):
    for j in range(i):
        dist = euclidean(physchem_descs[i], physchem_descs[j])
        if dist <= threshold:
            similarity = 1 - (dist - d_min) / (d_max - d_min)
            g.add_edge(smis[i], smis[j], weight=dist, label=f"{similarity:.2f}") 
            
            
# Define the edge style to display the edge labels (distances)

node_style_active = {"selector": "node[label = 'active']",
        "style": {
            "font-family": 'helvetica',
            "content": "data(name)",
            'width': 800,
            'height': 300,
            'shape': 'rectangle',
            'background-image': 'data(img)',
            "font-size": "60px",
            'border-opacity': 0,
            'border-width': 0.0,
            'background-fit': 'contain',
            "border-color": "green",  # Color for active compounds
            'background-color': 'green',
            "border-width": 3}
}
node_style_inactive = {       "selector": "node[label = 'inactive']",
        "style": {
            "font-family": 'helvetica',
            "content": "data(name)",
            'width': 800,
            'height': 300,
            'shape': 'circle',
            'background-image': 'data(img)',
            "font-size": "50px",
            'border-opacity': 0,
            'border-width': 0.0,
            'background-fit': 'contain',
            "border-color": "red",  # Color for inactive compounds
            'background-color': 'red',
            "border-width": 3
        }}

node_style_intermediate  = {  "selector": "node[label = 'intermediate']",
        "style": {
            "font-family": 'helvetica',
            "content": "data(name)",
            'width': 800,
            'height': 300,
            'shape': 'hexagon',
            'background-image': 'data(img)',
            "font-size": "50px",
            'border-opacity': 0,
            'border-width': 0.0,
            'background-fit': 'contain',
            "border-color": "blue",  # Color for inactive compounds
            'background-color': 'blue',
            "border-width": 3}
        }

edge_style = {
    "selector": "edge",
    "style": {
        "label": "data(label)",
        "font-size": "14px",
        "line-color": "grey",
        "width": 2,
    },
}            

Similarity_Graph= ipycytoscape.CytoscapeWidget()
Similarity_Graph.graph.add_graph_from_networkx(g)
styles = [node_style_active, node_style_inactive, node_style_intermediate, edge_style]
Similarity_Graph.set_style(styles)
Similarity_Graph.set_layout(name='dagre', padding=1,nodeSpacing=1,compactComponents= True,nodePlacement='SIMPLE',edgeSpacingFactor= 0.1)
display(Similarity_Graph) 
