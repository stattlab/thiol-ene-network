import sys, os, re
import numpy as np
import signac
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm
import gsd, gsd.hoomd 
import networkx as nx 
from collections import defaultdict
from matplotlib.pyplot import cm
import json


def connected_components(lists):
    R"""
    merges lists with common elements

    args: bonds from configuration (frame.bonds.group)

    returns: list of connected particles by id

    Useful for finding bonded particles in a configuration, tested for
    linear polymers with consecutive bonds (0-1-2-3-4-5, 6-7-8-9-10,..)
    and non consecutive ids ( 0-5-6-8-10, 1-4-3-2-9,...) but no other
    configuration yet. Works with ints as well as str.

    """
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

plt.rcParams["font.family"] = "Avenir"
fig, ax = plt.subplots(2,1,sharey=False)


project = signac.get_project()

color = iter(cm.rainbow(np.linspace(0, 1, 10)))

for job in project: 
    if job.isfile('polymerize.gsd'):
        # only look at our base cases for now...
        with open(job.fn('signac_statepoint.json')) as f:
            statepoint = json.load(f)
            if not (statepoint["crosslinker_percent"]==0.0 and statepoint["chain_side_reaction_probability"]==0):
                continue
        try:
            trajectory = gsd.hoomd.open(job.fn('polymerize.gsd'))
            frame = trajectory[0]
            dummy_id = len(frame.bonds.types)-1
            bonds = frame.bonds.group[frame.bonds.typeid!=dummy_id] # dummy type bond 
            nids, counts = np.unique(np.concatenate(bonds).flatten(),return_counts=True)
            all_strands = connected_components(bonds)
            n_molecules = 0
            for strand in all_strands:
                    n_molecules +=1
            N_0 = n_molecules
            c = next(color)
            conversion = []
            frame_number = []
            largest_molecule = []
            average_M = []
            for i,frame in enumerate(trajectory):
                bonds = frame.bonds.group[frame.bonds.typeid!=dummy_id]
                nids, counts = np.unique(np.concatenate(bonds).flatten(),return_counts=True)

                f_bonded = len(counts[counts==2])/frame.particles.N
                all_strands = connected_components(bonds)
                strand_lengths = []
                n_molecules = 0
                for strand in all_strands:
                    strand_lengths.append(len(strand))
                    n_molecules +=1
                extent_of_reaction = (N_0 - n_molecules)/N_0
                frame_number.append(i)
                conversion.append(extent_of_reaction)
                largest_molecule.append(np.max(strand_lengths))
                average_M.append(np.average(strand_lengths))

            conversion=np.array(conversion)
            ax[0].plot(conversion,2.0/(1.0-conversion),c='black')
            ax[0].plot(conversion,average_M, "o", c=c)
            
            ax[1].plot(frame_number,conversion,c=c)
            ax[0].set_ylim([1, 4.25])
            ax[1].set_ylim([0, 1])
            print(job.id, "done", len(trajectory))

        except:
            print("file exist but 0 size", job.id)
             

plt.show()