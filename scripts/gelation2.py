import sys, os, re, json, io, itertools
import numpy as np
import signac
#import matplotlib.pyplot as plt
#import matplotlib.cm
import gsd, gsd.hoomd 

import networkx as nx

def calc_weightAvgMw(G_list):
    num = 0
    denom = 0
    for G in G_list:
        M = len(G)
        num += M*M
        denom += M
    return num/denom
# G_list = [[1],[2,3],[3,3,3],[4,4,4,4],[5,5,5,5,5]]
# print(calc_weightAvgMw(G_list))

#---------------------------------------------------------------------------------------
# Analysis Functions
#---------------------------------------------------------------------------------------

#https://pubs.acs.org/doi/full/10.1021/acs.macromol.3c00831
# Calculates the reduced weight-average molecular weight of the system. This is the
# weight-avg Mw excluding the largest molecule.
def gelation_analysis_via_reducedMw(job_id):
    """
    Writes the gelation conversion based on the peak in reduced weight-average molecular
    weight to the job document as 'gelation_conversion_2_2'.
    """
    project = signac.get_project()
    job = project.open_job(id=job_id)
    traj = gsd.hoomd.open(job.fn('polymerize.gsd'))

    N0 = -1

    # iterate through the trajectory, and for each frame, calculate the reduced Mw
    reduced_Mws = []
    conversions = []
    for n,frame in enumerate(traj):
        bonds = frame.bonds.group
        G = nx.Graph()
        G.add_edges_from(bonds)
        Gcc = sorted(nx.connected_components(G), key=len, reverse=True)
        if n==0:
            N0 = len(Gcc)
        # make a list of molecules excluding the second largest one
        G_list = Gcc[1:]
        reduced_Mw = calc_weightAvgMw(G_list=G_list)
        reduced_Mws.append(reduced_Mw)
        conversions.append(1 - len(Gcc)/N0)
    
    #after iterating, find the conversion when the reduced Mw is at the peak
    max_idx = np.argmax(reduced_Mws)
    max_conv = conversions[max_idx]
    print("Conversion at reduced weight-averaged Mw peak:",max_conv)

    job.doc['gelation_conversion_2_2'] = (max_conv, max_idx)

    return max_conv


def gelation_analysis_via_2ndLargest(job_id):
    """
    Writes the gelation conversion based on the drop in size of the 2nd largest cluster
    to the job document as 'gelation_conversion_2'.
    """

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

    job.doc['gelation_conversion_2'] = (max_diff_conv, max_diff_idx)

    return max_diff_conv

def main(job_id):
    gelation_analysis_via_reducedMw(job_id)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 gelation2.py job_id")
        sys.exit(1)
    main(sys.argv[1])
