# handle image conversion
import networkx as nx
from networkx.readwrite import cytoscape_data
from rdkit import Chem
from rdkit.Chem import RDConfig
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs
from rdkit.Chem.Draw import rdMolDraw2D
import os
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
    #print(smi,nam)
    g.add_node(smi, img=smi2image(smi), name=nam, label=label)

# add edge if Tc >= threshold
mols = [Chem.MolFromSmiles(x) for x in smis]
fps = [AllChem.GetMorganFingerprintAsBitVect(x,3,2048) for x in mols]
for i in range(len(smis)):
    for j in range(i):
        Tc = DataStructs.TanimotoSimilarity(fps[i], fps[j])
        ### Pick the threshold default is 0.5
        if Tc >= 0.5:
            g.add_edge(smis[i], smis[j])

Similarity_Graph= ipycytoscape.CytoscapeWidget()
Similarity_Graph.graph.add_graph_from_networkx(g)
Similarity_Graph.set_style([
    {
        "selector": "node[label = 'active']",
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
            "border-width": 3
        }
    },
    {
        "selector": "node[label = 'inactive']",
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
        }
    },
    {
           "selector": "node[label = 'intermediate']",
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
            "border-width": 3
    }
    }
])


Similarity_Graph.set_layout(name='klay', padding=2,nodeSpacing=85)
display(Similarity_Graph) 
