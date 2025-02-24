import sys, os, re, json, io, itertools
import numpy as np
import signac
#import matplotlib.pyplot as plt
#import matplotlib.cm
import gsd, gsd.hoomd 

import networkx as nx


#---------------------------------------------------------------------------------------
# Analysis Functions
#---------------------------------------------------------------------------------------

def gelation_analysis_via_2ndLargest(job_id):
    #testing = True

    project = signac.get_project()
    job = project.open_job(id=job_id)
    traj = gsd.hoomd.open(job.fn('polymerize.gsd'))

    # get the number of molecules in the initial system
    N0 = -1

    # iterate through the trajectory, and for each frame, calculate the size of the 2nd largest cluster
    second_largest_cluster_sizes = []
    conversions = []
    for n,frame in enumerate(traj):
        bonds = frame.bonds.group
        G = nx.Graph()
        G.add_edges_from(bonds)
        Gcc = sorted(nx.connected_components(G), key=len, reverse=True)
        if n==0:
            N0 = len(Gcc)
        print(len(Gcc[1]))
        second_largest_cluster_sizes.append(len(Gcc[1]))
        conversions.append(1 - len(Gcc)/N0)
    
    #after iterating, find the conversion when the 2nd largest cluster size drops the most
    diff = np.diff(second_largest_cluster_sizes)
    max_diff_idx = np.argmax(diff)
    max_diff_conv = conversions[max_diff_idx]
    print("Conversion at 2nd largest cluster size drop:",max_diff_conv)
    '''
    if testing:
        ax, fig = plt.subplots()
        plt.plot(conversions, second_largest_cluster_sizes)
        plt.xlabel("Conversion")
        plt.ylabel("Size of 2nd Largest Cluster")
        plt.title("Conversion vs. Size of 2nd Largest Cluster")
        plt.show()
    '''
    job.doc['gelation_conversion_2'] = (max_diff_conv, max_diff_idx)

    return max_diff_conv

def main(job_id):
    gelation_analysis_via_2ndLargest(job_id)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 gelation2.py job_id")
        sys.exit(1)
    main(sys.argv[1])
