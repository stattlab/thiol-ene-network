import sys, os, re, json, io, itertools
import numpy as np
import signac
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
import matplotlib.pyplot as plt
import matplotlib.cm
=======
#import matplotlib.pyplot as plt
#import matplotlib.cm
>>>>>>> Stashed changes
=======
#import matplotlib.pyplot as plt
#import matplotlib.cm
>>>>>>> Stashed changes
=======
#import matplotlib.pyplot as plt
#import matplotlib.cm
>>>>>>> Stashed changes
import gsd, gsd.hoomd 

import networkx as nx


#---------------------------------------------------------------------------------------
# Analysis Functions
#---------------------------------------------------------------------------------------

def gelation_analysis_via_2ndLargest(job_id):
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
    testing = True
=======
    #testing = True
>>>>>>> Stashed changes
=======
    #testing = True
>>>>>>> Stashed changes
=======
    #testing = True
>>>>>>> Stashed changes

    project = signac.get_project()
    job = project.open_job(id=job_id)
    traj = gsd.hoomd.open(job.fn('polymerize.gsd'))

    # get the number of molecules in the initial system
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
    N0 = job.sp["N_monomers"]
=======
    N0 = -1
>>>>>>> Stashed changes
=======
    N0 = -1
>>>>>>> Stashed changes
=======
    N0 = -1
>>>>>>> Stashed changes

    # iterate through the trajectory, and for each frame, calculate the size of the 2nd largest cluster
    second_largest_cluster_sizes = []
    conversions = []
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
    for frame in traj:
=======
    for n,frame in enumerate(traj):
>>>>>>> Stashed changes
=======
    for n,frame in enumerate(traj):
>>>>>>> Stashed changes
=======
    for n,frame in enumerate(traj):
>>>>>>> Stashed changes
        bonds = frame.bonds.group
        G = nx.Graph()
        G.add_edges_from(bonds)
        Gcc = sorted(nx.connected_components(G), key=len, reverse=True)
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
        print(len(Gcc[1].nodes()))
        second_largest_cluster_sizes.append(len(Gcc[1].nodes()))
=======
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
        if n==0:
            N0 = len(Gcc)
        print(len(Gcc[1]))
        second_largest_cluster_sizes.append(len(Gcc[1]))
<<<<<<< Updated upstream
<<<<<<< Updated upstream
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
        conversions.append(1 - len(Gcc)/N0)
    
    #after iterating, find the conversion when the 2nd largest cluster size drops the most
    diff = np.diff(second_largest_cluster_sizes)
    max_diff_idx = np.argmax(diff)
    max_diff_conv = conversions[max_diff_idx]
    print("Conversion at 2nd largest cluster size drop:",max_diff_conv)
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream

=======
    '''
>>>>>>> Stashed changes
=======
    '''
>>>>>>> Stashed changes
=======
    '''
>>>>>>> Stashed changes
    if testing:
        ax, fig = plt.subplots()
        plt.plot(conversions, second_largest_cluster_sizes)
        plt.xlabel("Conversion")
        plt.ylabel("Size of 2nd Largest Cluster")
        plt.title("Conversion vs. Size of 2nd Largest Cluster")
        plt.show()
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream

    job.doc['gelation_conversion_2'] = max_diff_conv
=======
    '''
    job.doc['gelation_conversion_2'] = (max_diff_conv, max_diff_idx)
>>>>>>> Stashed changes
=======
    '''
    job.doc['gelation_conversion_2'] = (max_diff_conv, max_diff_idx)
>>>>>>> Stashed changes
=======
    '''
    job.doc['gelation_conversion_2'] = (max_diff_conv, max_diff_idx)
>>>>>>> Stashed changes

    return max_diff_conv

def main(job_id):
    gelation_analysis_via_2ndLargest(job_id)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 gelation2.py job_id")
        sys.exit(1)
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
    main(sys.argv[1])
=======
    main(sys.argv[1])
>>>>>>> Stashed changes
=======
    main(sys.argv[1])
>>>>>>> Stashed changes
=======
    main(sys.argv[1])
>>>>>>> Stashed changes
