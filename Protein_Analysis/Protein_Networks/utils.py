import numpy as np
from networkx import from_numpy_matrix, to_numpy_matrix
from urllib.request import urlretrieve
from scipy.spatial.distance import euclidean
from operator import itemgetter
#centrality measures
from networkx.algorithms.centrality import degree_centrality, closeness_centrality, betweenness_centrality, betweenness_centrality_subset
import datetime
import os

def readPDBFile(pdbFilePath):

    datetime_object = datetime.datetime.now()
    print('Start Reading PDB'+'\n')
    print(datetime_object)
    atoms = []
    with open(pdbFilePath) as pdbfile:
        for line in pdbfile:
            if line[:4] == 'ATOM':
                #split the line
                splitted_line = [line[:6], line[6:11], line[12:16], line[17:20], line[21], line[22:26], line[30:38], line[38:46], line[46:54], line[56:61], line[62:66]]
                atoms.append(splitted_line)
    datetime_object = datetime.datetime.now()
    print('End Reading PDB'+'\n')
    print(datetime_object)
    return np.array(atoms)

def getResidueCoordinates(atoms):

    coordinates = []
    residues_list = ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS',
                   'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL',
                   'TRP', 'TYR']

    dropped = []

    last_residue_num = 0

    for i, atom in enumerate(atoms):

        residue_name = atom[3]
        residue_chain = atom[4]
        residue_num = int (atom[5])

        if residue_name in residues_list:

            if (residue_num != last_residue_num):

                if (atom[2].replace(" ","")=="CA"):
                    cord_C = atom[6:9]
                    coordinates.append([residue_name + str (residue_num) + " " + residue_chain, cord_C])
                    last_residue_num = residue_num

        else:
          dropped.append(residue_name)

    return np.array(coordinates)

def associateResidueName(coordinates):

    dict_residue_name = dict()
    for i in range(coordinates.shape[0]):
        dict_residue_name[str (i)] = coordinates[i, 0]
    return dict_residue_name

def getResiduesSequence(pbdFilePath):

    seq_res = []
    residues_list = ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR']

    with open(pbdFilePath) as pdbfile:
        for line in pdbfile:
            if line[:6]=='SEQRES':
                splitted_line = [line[19:22], line[23:26], line[27:30], line[31:34], line[35:38],
                                line[39:42], line[43:46], line[47:50], line[51:54], line[55:58],
                                line[59:62], line[63:66], line[67:70]]
                for residue in splitted_line:
                    if residue in residues_list:
                        seq_res.append(residue)

    return np.array(seq_res)

def printProgressBar (iteration, total):

    length = 100
    fill = '█'
    printEnd = '\r'
    percent = ("{0:." + str(2) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print(f"\r |{bar}| Current progress: {percent}%", end = printEnd)

    if percent == 100:
        print()

def adjacent_matrix(output_path, coordinates, p, min_=4, max_=8):
    from sklearn.metrics import pairwise_distances
    import time

    start = time.time()
    n = coordinates.shape[0]
    cords = np.hstack(coordinates[:,1]).reshape(n,3)
    adj = np.zeros((n,n))
    #d = np.zeros((n,n), dtype=float)
    d = pairwise_distances(X = cords, n_jobs = -1)

    # if comp_adj_fr is not None:
    #     pb = ttk.Progressbar(comp_adj_fr, orient="horizontal", mode="determinate", length=100)
    #     pb.pack()
    #     pb["value"] = 0
    #     label = tk.Label(comp_adj_fr, text="Current progress {}%".format(pb["value"]))
    #     label.pack()
    #     window.update()
    # else:
    value = 0
    printProgressBar(value, n)

    for i in range(n):
        for j in range(i):
            if ((d[i][j]>min_) and (d[i][j]<max_)):
                adj[i][j] = 1
                adj[j][i] = 1

        # if comp_adj_fr is not None:
        #     pb["value"] = round(((i + 1) / n) * 100, 2)
        #     label['text'] = "Current progress {}%".format(pb["value"])
        #     pb.pack()
        #     label.pack()
        #     window.update()
        # else:
        #     printProgressBar(i + 1, n)

    # if comp_adj_fr is not None:
    #     pb["value"] = round(((i + 1) / n) * 100, 2)
    #     label['text'] = "Current progress {}%".format(pb["value"])
    #     pb.pack()
    #     label.pack()
    #     window.update()
    # else:
    #    # printProgressBar(i + 1, n)

    end = time.time()
    print("Time for parallel PCN computation of protein {}: {} s".format(p, (end-start)))

    if not os.path.exists("{}Adj".format(output_path)):
        os.makedirs("{}Adj".format(output_path))
    np.savetxt("{}Adj{}_adj_mat_{}_{}.txt".format(output_path, p, min_, max_), adj, fmt='%.2f')
    print("saved adj matrix")

    return adj

def degree_centrality(G, res_names, n=10):
    dc = degree_centrality(G)
    dc = {int (float (k)):v for k,v in dc.items()}
    dict_node_centrality = dict ()

    for i, cent in dc.items():

        dict_node_centrality[res_names[i]] = cent

    sorted_dc = sorted(dict_node_centrality.items(), key=itemgetter(1), reverse=True)
    print("Top {} nodes by degree centrality".format(n))
    for d in sorted_dc[:n]:
        print(d)

    return dict_node_centrality

def closeness(G, res_names, n=10):
   
    cc = closeness_centrality(G)
    cc = {int (float (k)):v for k,v in cc.items()}
    dict_node_centrality = dict ()

    for i, cent in cc.items():

        dict_node_centrality[res_names[i]] = cent

    sorted_cc = sorted(dict_node_centrality.items(), key=itemgetter(1), reverse=True)
    print("Top {} nodes by closeness_centrality".format(n))
    for d in sorted_cc[:n]:
        print(d)

    return dict_node_centrality


def betweenness(G, res_names, n=10):

    bc = betweenness_centrality(G)
    #bc= betweenness_centrality_parallel(G)
    bc = {int (float (k)):v for k,v in bc.items()}
    dict_node_centrality = dict ()

    for i, cent in bc.items():

        dict_node_centrality[res_names[i]] = cent

    sorted_bc = sorted(dict_node_centrality.items(), key=itemgetter(1), reverse=True)
    print("Top {} nodes by betweenness centrality".format(n))
    for d in sorted_bc[:n]:
        print(d)

    return dict_node_centrality
