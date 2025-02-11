import sys, os, re, json, io, itertools
import numpy as np
import pandas as pd
import signac
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm
from matplotlib.pyplot import cm
from matplotlib.colors import rgb2hex
import gsd, gsd.hoomd 
from collections import defaultdict
from collections import OrderedDict
from collections import Counter

import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import seaborn as sns
import plotly.io as pio   

import networkx as nx

pio.kaleido.scope.mathjax = None

project = signac.get_project()
direc = project.fn('') + 'workspace/'

"""
frame.bonds.group gives a np.array that is of form [[a,b],[c,d],[e,f],[a,g],...] 
where a is bonded to b, c is bonded to d, etc.

atom_bonds finds all atoms to which a given atom 'a' is bonded to
args: frame.bonds.group

returns: dictionary of atoms and the atoms they are bonded to
{atom_no: [bonded1, bonded2, ...], ...}
"""
def atom_bonds(bonds):
    # create a default dictionary with an element returning a set if nothing in it
    neighbors = defaultdict(list) 
    # add every instance of an atom in a bond to the dictionary
    for each in bonds:
        for item in each:
            neighbors[item].extend(each)
            # remove the atom itself from the list of atoms bonded to it
            neighbors[item].remove(item)
    return neighbors


def directly_bonded(beads, bonds):
    direct_bonds = []
    for each in bonds:
        for n,bead in enumerate(each):
            if bead in beads:
                direct_bonds.append(each[1-n])
                break
    return direct_bonds



for job_id in os.listdir(direc):
    
    jdir = direc + job_id
    input_file = direc + job_id + "/polymerize.gsd"
    if not os.path.isfile(input_file):
        continue
    # open the polymerized gsd file's trajectory (list of frames)
    trajectory = gsd.hoomd.open(input_file)

    frame = trajectory[-1]
    dummy_id = len(frame.bonds.types)-1
    bonds = frame.bonds.group[frame.bonds.typeid!=dummy_id]


    Ntotal = frame.particles.N

    # make a graph - each bond is an edge
    G = nx.Graph()
    G.add_edges_from(bonds)
    # find all disconnected/connected sub-networks
    # sort the list of networks by length
    Gcc = sorted(nx.connected_components(G), key=len, reverse=True)

    '''
    #print("number of disconnected clusters in network",len(Gcc))
    sizes = [len(n) for n in Gcc ]
    #print("sizes of clusters in network",sizes)
    
    # calculate gel fraction
    sizes = sorted(Counter(sizes).items(), key=lambda item: item[0], reverse=True)
    with open(jdir + "/cluster_sizes.txt", "w") as f:
        f.write("cluster_size count\n")
        for x in sizes:
            line = str(x[0]) + " " + str(x[1]) + "\n"
            f.write(line)
    '''
        
    # strand lengths (in entrie system. If you only want the ones in the gel, you need to look at Gcc[0] only)
    # bonds in biggest cluster:
    new_bonds = []
    for each in bonds:
        for i in each:
            if i in Gcc[0]:
                new_bonds.append(each)
                continue
    # make cluster of biggest graph
    G2 = nx.Graph()
    G2.add_edges_from(new_bonds)

    '''
    Gnew = G.copy()
    # remove everyone that has 3 or more bonds on it, only leaving linear strands
    crosslink_beads = [x for  x in G2.nodes() if G2.degree(x) >= 3]


    ## find beads that are directly attached to crosslinks
    print("finding beads directly bonded to crosslinks")
    beads_bonded_to_crosslinks = np.unique(directly_bonded(crosslink_beads, new_bonds))

    
    # all_bonds = atom_bonds(new_bonds)
    # print("beep")
    # print(all_bonds)    
    # print("beep")
    # beads_bonded_to_crosslinks2 = [all_bonds[atom] for atom in all_bonds.keys() if len(all_bonds[atom])>2]
    # beads_bonded_to_crosslinks2 = np.unique([x for xs in beads_bonded_to_crosslinks2 for x in xs])

    # print(len(set(beads_bonded_to_crosslinks2).intersection(set(beads_bonded_to_crosslinks))) == len(beads_bonded_to_crosslinks2))
    # print(len(set(beads_bonded_to_crosslinks2).intersection(set(beads_bonded_to_crosslinks))))
    # print(len(beads_bonded_to_crosslinks))
    # print(len(beads_bonded_to_crosslinks2))

    for x in crosslink_beads:
        Gnew.remove_node(x)

    # strands
    Gcc_new = sorted(nx.connected_components(Gnew), key=len, reverse=True)
    print("number of strands network",len(Gcc_new))
    sizes = [len(n) for n in Gcc_new]
    print("length of strands in network")
    sizes_of_strands,count = np.unique(sizes,return_counts=True)
    print(sizes_of_strands,count)

    strand_count = list(zip(sizes_of_strands,count))
    print(strand_count)
    with open(jdir + "/strand_sizes.txt", "w") as f:
        f.write("strand_size count\n")
        for x in strand_count:
            line = str(x[0]) + " " + str(x[1]) + "\n"
            f.write(line)
    
    # dangling ends
    print("mapping dangling ends")
    dangling_ends = []
    for i in Gcc_new:
        intersect = i.intersection(set(beads_bonded_to_crosslinks))
        if len(intersect) < 2 and len(intersect) != len(i):
            dangling_ends.append(len(i))
    
    sizes_of_dangles,count = np.unique(dangling_ends,return_counts=True)
    dangle_count = list(zip(sizes_of_dangles,count))
    print(dangle_count)
    with open(jdir + "/dangle_sizes.txt", "w") as f:
        f.write("dangle_size count\n")
        for x in dangle_count:
            line = str(x[0]) + " " + str(x[1]) + "\n"
            f.write(line)
    '''
    
    print("calculating loops")
    loops = list(nx.simple_cycles(G2,length_bound=50))
    lengths_loops = [len(l) for l in loops]
    sizes_of_loops,count = np.unique(lengths_loops,return_counts=True)
    loop_count = list(zip(sizes_of_loops,count))
    print(loop_count)
    with open(jdir + "/loop_counts.txt", "w") as f:
        f.write("loop_size count\n")
        for x in loop_count:
            line = str(x[0]) + " " + str(x[1]) + "\n"
            f.write(line)
