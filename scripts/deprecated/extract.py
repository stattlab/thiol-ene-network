import sys, os, re
import numpy as np
import signac
import json
import gsd, gsd.hoomd 
import networkx as nx 
from collections import defaultdict



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

"""
merges lists with common elements

args: bonds from configuration (frame.bonds.group)

returns: list of connected particles by id

Useful for finding bonded particles in a configuration, tested for
linear polymers with consecutive bonds (0-1-2-3-4-5, 6-7-8-9-10,..)
and non consecutive ids ( 0-5-6-8-10, 1-4-3-2-9,...) but no other
configuration yet. Works with ints as well as str.

"""
def connected_components(lists):
    
    neighbors = defaultdict(set)
    seen = set()
    for each in lists:
        for item in each:
            neighbors[item].update(each)
    def component(node, neighbors=neighbors, seen=seen, see=seen.add):
        nodes = set([node])
        next_node = nodes.pop
        while nodes:
            node = next_node()
            see(node)
            nodes |= neighbors[node] - seen
            yield node
    for node in neighbors:
        if node not in seen:
            yield sorted(component(node))

### a---b
###  \ /
###   c
"""
Counts and identifies cyclopropyl groups in the system
args: bonds from configuration (frame.bonds.group)
"""
def cyclopropyl_groups(a_bonds):
    cyclopropyl_count = 0
    # get keys for atoms with more than one bond
    keys = [i for i in a_bonds.keys() if len(a_bonds[i])>1]
    cyclo_keys = []
    cyclo_keys_flat = []
    # a
    for a in keys:
        if a in cyclo_keys_flat:
            continue
        # b
        for b in a_bonds[a]:
            if b in cyclo_keys_flat:
                continue
            # c
            for c in a_bonds[b]:
                if c in cyclo_keys_flat:
                    continue
                if a in a_bonds[c]:
                    cyclo_keys = cyclo_keys + [[a,b,c]]
                    cyclopropyl_count += 1
                    cyclo_keys_flat = [x for xs in cyclo_keys for x in xs]
    return(cyclopropyl_count, cyclo_keys)


''' find amount of atoms with x number of bonds'''
def count_atom_bonds(all_bonds):
    # count the atoms with (1,2,3,4, or 5) bonds
    bond_counts = defaultdict(int, {k:0 for k in (1,2,3,4) })
    for atom, neighbors in all_bonds.items():
        bond_counts[int(len(neighbors))] += 1 
    return bond_counts

''' GET ALL BONDS in the biggest molecules
find amount of crosslinks that occur from thiol ene reaction and those that occur from chain growht
crosslink defined as molecule bonded to at least 3 different starting reagents
'''
def crosslink_types(biggest_molec_all_bonds):
    # go through all atoms and their bonds to see what they are bonded to
    crosslinks = {"thiol-ene":0, "chain-growth":0}
    for atom, neighbors in biggest_molec_all_bonds.items():
        ## CHECKING TYPE OF CROSSLINK
        # THIOL-ENE - 4 but 3 of the 4 arms must be reacted
        if len(neighbors) == 4:
            arms_reacted = 0
            for n in neighbors:
                if len(all_bonds[n]) > 0:
                    arms_reacted += 1
            if arms_reacted > 3:
                crosslinks["thiol-ene"] += 1
        # CHAIN-GROWTH if it has 3 bonds
        if len(neighbors) == 3:
            crosslinks["chain-growth"] += 1
    return crosslinks


# count the types of each unreacted atom:
def unreacted_atoms(type_map):
    non_reacted_atoms = [types[j] for j in type_map if j<4]
    unreact_hist = {'Thiol':0,'Ene_C':0,'RSulfur':0,'RCarbon':0}
    for j in non_reacted_atoms:
        unreact_hist[j] = unreact_hist.get(j, 0) + 1
    return(unreact_hist)


def parse_final_frame(input_file):
    # open the polymerized gsd file's trajectory (list of frames)
    trajectory = gsd.hoomd.open(input_file)
    # get the dummy_bond_id and number of bonds from the initial frame
    initial_frame = trajectory[0]
    dummy_id = len(frame.bonds.types)-1
    bonds = frame.bonds.group[frame.bonds.typeid!=dummy_id] # dummy type bond 
    # find all strands in the box
    all_strands = connected_components(bonds)
    # count the number of molecules
    n_molecules = 0
    for strand in all_strands:
            if cg_only and len(strand) == 2:
                continue
            n_molecules +=1
    N_0 = n_molecules

    # conversion, largest molecule size and average molecule size through the trajectory

    
    ## initialize information
    # look at the final frame
    frame = trajectory[-1]
    # get non_dummy_bonds
    bonds = frame.bonds.group[frame.bonds.typeid!=dummy_id]
    Ntotal = frame.particles.N
    # atoms and what they are connected to
    all_bonds = atom_bonds(bonds)
    # all strings of atoms
    all_strands = connected_components(bonds)
    # the types of atoms
    type_map = frame.particles.typeid


    ## analyze
    # get the count of atoms with x number bonds
    bond_counts = count_atom_bonds(all_bonds)
    # counts of reactions produced each crosslink
    CLtype_counts = crosslink_types(biggest_molec_all_bonds)
    # counts of unreacted atoms
    unreacted_counts = unreacted_atoms(type_map)


    # network analysis
    # make a graph - each bond is an edge
    G = nx.Graph()
    G.add_edges_from(bonds)
    G.add_edges_from(bonds)

    # find all disconnected/connected sub-networks
    # sort the list of networks by length
    Gcc = sorted(nx.connected_components(G), key=len, reverse=True)
    print("number of disconnected clusters in network",len(Gcc))
    sizes = [len(n) for n in Gcc ]
    print("sizes of clusters in network",sizes)
    

direc = ""


