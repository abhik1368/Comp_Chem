from __future__ import print_function
import numpy as np
from pymol import cmd
from pymol.querying import get_color_indices 
import os 

add_slash_to_path ='\\'
def get_colors(selection='', quiet=1):
    
    pymol_color_list = []
    
    for tuplepair in get_color_indices(selection):
        pymol_color_list.append(tuplepair[0])
    
    #pymol_color_list.sort()
    if not int(quiet): print(pymol_color_list)
    pymol_color_list.remove('black')
    pymol_color_list.remove('white')
    pymol_color_list.remove('dash')
    return pymol_color_list

cmd.extend('get_colors',get_colors)
cmd.auto_arg[0]['get_colors']=[lambda: cmd.Shortcut(['""','all']), 'selection=', ',']
cmd.auto_arg[1]['get_colors']=[lambda: cmd.Shortcut(['0']), 'quiet=', '']

def pymol_centralities(output_path, centralities, protein_path, algorithm_name):
    
    cmd.do("delete {}".format("all"))
    cmd.do("load {}".format(protein_path))
    cmd.do("set specular, off")
    protein = os.path.basename(protein_path)
    protein_name = os.path.splitext(protein)[0]
    
    colors = get_colors()
    
    cmd.do("remove hetatm")
    
        
    n = len(list(centralities.keys()))
    for count, (residue, cent) in enumerate(centralities.items()):
        residue_n, residue_chain = residue.split(" ")
        residue_name = residue_n[:3]
        residue_num = residue_n[3:]
        line="alter (resi "+ str(residue_num) + " and chain "+ residue_chain + "), b = "+ str(cent)
        cmd.do(line)           
    cmd.do("spectrum b, rainbow")
    cmd.do("ramp_new colorbar, none, [{}, {}], rainbow".format(min(centralities.values()), max(centralities.values())))
    
    if (not os.path.exists("{}Centralities{}{}{}Sessions".format(output_path, add_slash_to_path, algorithm_name, add_slash_to_path))):
        os.makedirs("{}Centralities{}{}{}Sessions".format(output_path, add_slash_to_path, algorithm_name, add_slash_to_path))
        
    cmd.do("save {}Centralities{}{}{}Sessions{}{}_{}_session.pse".format(output_path, add_slash_to_path, algorithm_name, add_slash_to_path, add_slash_to_path, protein_name, algorithm_name))

