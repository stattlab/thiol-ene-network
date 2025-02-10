import sys, os, re
import numpy as np
import signac
import json
import gsd, gsd.hoomd 
import networkx as nx 
from collections import defaultdict

from extract import connected_components

'''
VIA GSD.HOOMD
for every frame of the polymerization trajectory
by going through each string of connected atoms, calculate conversion, 
the largest molecule by bead size, and the average molecule bead size 
'''
def conversion_molecule_sizes(N_0, trajectory):
    traj_conversion = []
    traj_molecule_sizes = []
    for i,frame in enumerate(trajectory):
        bonds = frame.bonds.group[frame.bonds.typeid!=dummy_id]
        connections = connected_components(bonds)
        strand_lengths = []
        n_molecules = 0
        for strand in connections:
            strand_lengths.append(len(strand))
            n_molecules +=1

        # the extent of the reaction
        extent_of_reaction = (N_0 - n_molecules)/N_0
        traj_conversion.append((i,extent_of_reaction))

        # get the largest molecule size and the average molecule size
        traj_molecule_sizes.append((np.max(strand_lengths), np.average(strand_lengths)))
        print(extent_of_reaction, np.max(strand_lengths))
        if i == len(trajectory) - 1:
            # count the number of each connection with size N
            molec_size_hist = dict()
            for j in strand_lengths:
                molec_size_hist[j] = molec_size_hist.get(j, 0) + 1
    return traj_conversion, traj_molecule_sizes, molec_size_hist

# direc = "/Users/rithwikghanta/Documents/Documents/research/network_polymerization/workspace/"
project = signac.get_project()
direc = project.fn('') + 'workspace/'
for job_id in os.listdir(direc):
    print(job_id)
    jdir = direc + job_id
    input_file = direc + job_id + "/polymerize.gsd"
    if not os.path.isfile(input_file):
        continue
    
    # open the polymerized gsd file's trajectory (list of frames)
    trajectory = gsd.hoomd.open(input_file)
    # get the dummy_bond_id and number of bonds from the initial frame
    initial_frame = trajectory[0]
    dummy_id = len(initial_frame.bonds.types)-1
    bonds = initial_frame.bonds.group[initial_frame.bonds.typeid!=dummy_id] # dummy type bond 
    # find all strands in the box
    all_strands = connected_components(bonds)
    # count the number of molecules
    n_molecules = 0
    for strand in all_strands:
            n_molecules +=1
    N_0 = n_molecules

    t_c, t_m_s, m_s_h_final = conversion_molecule_sizes(N_0, trajectory)

    # write trajectory conversion into file
    with open(jdir + "/trajectory_conversion.txt", "w") as f:
        f.write("frame conversion\n")
        for x in t_c:
            line = str(x[0]) + " " + str(x[1]) + "\n"
            f.write(line)

    # write trajectory largest and avg molecule sizes into file
    with open(jdir + "/trajectory_molecule_size.txt", "w") as f:
        f.write("largest average\n")
        for x in t_m_s:
            line = str(x[0]) + " " + str(x[1]) + "\n"
            f.write(line)

    # write final frame molecule size counts into file
    with open(jdir + '/molecule_size_histogram.json', "w") as f:
        json.dump(m_s_h_final, f)