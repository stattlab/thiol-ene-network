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
from datetime import datetime

import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import seaborn as sns
import plotly.io as pio   

import networkx as nx

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

def classify_loop_type(cycle,bonds):
    n_multifunctional = 0
    for bead in cycle:
        n_bonds = 0
        for bond in bonds:
            if bead in bond:
                n_bonds += 1
        if n_bonds > 2:
            n_multifunctional += 1
    if n_multifunctional == 0:
        # return "primary"
        raise ValueError("Complete loop found. Are you running on the largest cluster?")
    elif n_multifunctional == 1:
        return "primary"
    elif n_multifunctional == 2:
        return "secondary"
    elif n_multifunctional == 3:
        return "tertiary"
    else:
        return "quaternary+"

def snap_molecule_indices(snap):
    """Find molecule index for each particle.

    Given a snapshot from a trajectory, compute clusters of bonded molecules
    and return an array of the molecule index of each particle.

    Parameters
    ----------
    snap : gsd.hoomd.Snapshot
        Trajectory snapshot.

    Returns
    -------
    numpy array (N_particles,)

    """
    system = freud.AABBQuery.from_system(snap)
    num_query_points = num_points = snap.particles.N
    query_point_indices = snap.bonds.group[:, 0]
    point_indices = snap.bonds.group[:, 1]
    # print(np.shape(query_point_indices))
    # print(np.shape(point_indices))
    # This may fix the error of the indices not being sorted, but uncertain
    # sorted_indices = np.argsort(query_point_indices)
    # query_point_indices = query_point_indices[sorted_indices]
    # point_indices = point_indices[sorted_indices]
    vectors = system.box.wrap(
        system.points[query_point_indices] - system.points[point_indices]
    )
    nlist = freud.NeighborList.from_arrays(
        num_query_points, num_points, query_point_indices, point_indices, vectors
    )
    cluster = freud.cluster.Cluster()
    cluster.compute(system=system, neighbors=nlist)
    return cluster.cluster_idx

def intermolecular_rdf(
    gsd_traj,
    A_name=None,
    B_name=None,
    r_max=None,
    r_min=0,
    bins=100,
    exclude_bonded=True,
):
    """Compute intermolecular RDF from a GSD file.

    This function calculates the radial distribution function given a GSD frame
    and the names of the particle types. By default it will calculate the RDF
    for all particles.

    It is assumed that the bonding, number of particles, and simulation box do
    not change during the simulation.

    Parameters
    ----------
    gsd_traj : gsd.hoomd.Trajectory
        The GSD trajectory.
    A_name, B_name : str
        Name(s) of particles between which to calculate the RDF (found in
        gsd.hoomd.Snapshot.particles.types)
    r_max : float
        Maximum radius of RDF. If None, half of the maximum box size is used.
        (Default value = None)
    r_min : float
        Minimum radius of RDF. (Default value = 0)
    bins : int
        Number of bins to use when calculating the RDF. (Default value = 100)
    exclude_bonded : bool
        Whether to remove particles in same molecule from the neighbor list.
        (Default value = True)

    Returns
    -------
    freud.density.RDF
    """
    import freud
    
    snap = gsd_traj[0]

    if r_max is None:
        # Use a value just less than half the maximum box length.
        r_max = np.nextafter(
            np.max(snap.configuration.box[:3]) * 0.5, 0, dtype=np.float32
        )

    rdf = freud.density.RDF(bins=bins, r_max=r_max, r_min=r_min)

    if A_name == None and B_name == None:
        type_A = [True]*snap.particles.N
        type_B = [True]*snap.particles.N
    else:
        type_A = snap.particles.typeid == snap.particles.types.index(A_name)
        type_B = snap.particles.typeid == snap.particles.types.index(B_name)

    if exclude_bonded:
        molecules = snap_molecule_indices(snap)
        molecules_A = molecules[type_A]
        molecules_B = molecules[type_B]

    for snap in gsd_traj:

        # operate on the frame
        A_pos = snap.particles.position[type_A]
        if A_name == B_name:
            B_pos = A_pos
            exclude_ii = True
        else:
            B_pos = snap.particles.position[type_B]
            exclude_ii = False

        print(A_pos)
        print(np.shape(A_pos))

        box = snap.configuration.box
        system = (box, A_pos)
        aq = freud.locality.AABBQuery.from_system(system)
        nlist = aq.query(
            B_pos, {"r_max": r_max, "exclude_ii": exclude_ii}
        ).toNeighborList()

        if exclude_bonded:
            pre_filter = len(nlist)
            indices_A = molecules_A[nlist.point_indices]
            indices_B = molecules_B[nlist.query_point_indices]
            nlist.filter(indices_A != indices_B)
            post_filter = len(nlist)

        rdf.compute(aq, neighbors=nlist, reset=False)
    normalization = post_filter / pre_filter if exclude_bonded else 1
    return rdf, normalization

def reduce_network(in_graph):
    dangling_ends = []
    out_graph = in_graph.copy()
    for i in in_graph.nodes():
        if in_graph.degree(i) == 1:
            out_graph.remove_node(i)
            dangling_ends.append(i)
    return out_graph, dangling_ends

def remove_dangling_ends(in_graph):
    dangling_ends = []
    new_graph = nx.Graph()
    last_graph = in_graph.copy()
    for i in range(1000):
        new_graph, new_dangling_ends = reduce_network(last_graph)
        dangling_ends.extend(new_dangling_ends)
        # print("last graph:",last_graph.edges())
        # print("new graph:",new_graph.edges())
        if new_graph.edges() == last_graph.edges():
            break
        else:
            last_graph = new_graph.copy()
    return new_graph, dangling_ends
    

#---------------------------------------------------------------------------------------
# Analysis Functions
#---------------------------------------------------------------------------------------

#Dangling end lengths do not include the beads within the gel
#
def defect_analysis(job_id):
    testing = False
    
    pio.kaleido.scope.mathjax = None

    project = signac.get_project()
    direc = project.fn('') + 'workspace/'

    # for job_id in os.listdir(direc):
    jdir = direc + job_id
    input_file = direc + job_id + "/polymerize_trunc.gsd"
    if not os.path.isfile(input_file):
        raise FileNotFoundError("File not found")
        exit()
    # open the polymerized gsd file's trajectory (list of frames)
    trajectory = gsd.hoomd.open(input_file)

    frame = trajectory[-1]
    dummy_id = len(frame.bonds.types)-1
    bonds = frame.bonds.group[frame.bonds.typeid!=dummy_id]


    Ntotal = frame.particles.N

    if testing:
        jdir = "./test_network_properties"
        bonds = [[0,1],[1,2],[2,3],[3,4],[4,1],[0,5],[5,6],[6,7],[5,8],[8,9],
                    [2,10],[2,11],[10,12],[12,11],[12,13],[12,14],[13,15],[15,14]]
        G = nx.Graph()
        G.add_edges_from(bonds)
        G2 = nx.Graph()
        G2.add_edges_from(bonds)
        frame = gsd.hoomd.Frame()
        frame.particles.position = np.array([[0,0,0],[1,0,0],[2,0,0],[2,1,0],[1,1,0],[-1,0,0],[-2,0,0],[-3,0,0],[-1,1,0],[-1,2,0],
                                             [3,0,0],[2,1,0],[3,1,0],[4,1,0],[3,2,0],[4,2,0]])
        frame.particles.velocity = np.zeros((len(frame.particles.position),3))
        frame.particles.N = len(frame.particles.position)
        frame.bonds.N = len(bonds)
        frame.bonds.group = np.array(bonds)
        frame.bonds.typeid = np.zeros(len(bonds))
        frame.bonds.types = ['bond']
        new_bonds = bonds

    else:
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
            
        # strand lengths (in entire system. If you only want the ones in the gel, you need to look at Gcc[0] only)
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
    
    print("calculating dangling ends")
    dangling_end_strands = []
    dangling_end_lengths = []

    reduced_graph, dangling_ends = remove_dangling_ends(G2)
    dangling_end_graph = G2.copy()
    for x in G2.nodes():
        if x not in dangling_ends:
            dangling_end_graph.remove_node(x)
    print("dangling ends:",dangling_end_graph.nodes())
    print("dangling end edges:",dangling_end_graph.edges())

    for dangling_strand in nx.connected_components(dangling_end_graph):
        dangling_strand = list(dangling_strand)
        print(dangling_strand)
        dangling_end_strands.append(dangling_strand)
        dangling_end_lengths.append(len(dangling_strand))
    unique_dangling_end_lengths,count = np.unique(dangling_end_lengths,return_counts=True)
    dangling_end_data = list(zip(unique_dangling_end_lengths,count))
    print(dangling_end_data)

    '''
    # This was the old way of calculating dangling ends, removed to allow for branching 
    #     dangling ends
    G2_strands = G2.copy()
    # remove everyone that has 3 or more bonds on it, only leaving linear strands
    crosslink_beads = [x for  x in G2.nodes() if G2.degree(x) >= 3]
    print("crosslink_beads:",crosslink_beads)
    for x in crosslink_beads:
        G2_strands.remove_node(x)
    
    for strand in nx.connected_components(G2_strands):
        strand = list(strand)
        # intersect = set(strand).intersection(set(crosslink_beads))
        # print("intersect:",intersect)
        # Checked and there is no intersection!
        print("strand:",strand)
        if G.degree(strand[0]) == 1 or G.degree(strand[-1]) == 1:
            dangling_end_strands.append(strand)
            dangling_end_lengths.append(len(strand))
    unique_dangling_end_lengths,count = np.unique(dangling_end_lengths,return_counts=True)
    dangling_end_data = list(zip(unique_dangling_end_lengths,count))
    '''

    # if (testing):
    #     out_frame = frame
    #     out_frame.particles.velocity[:] = -1
    #     for i,strand in enumerate(nx.connected_components(G2_strands)):
    #         strand = list(strand)
    #         for bead in strand:
    #             out_frame.particles.velocity[bead] = [i,0,0]
    #     output_gsd = gsd.hoomd.open(jdir + f'/test_{os.path.basename(__file__)}'.replace('.py','.gsd'), 'w')
    #     output_gsd.append(out_frame)
    #     exit()

    print("calculating loops")
    start = datetime.now()
    print('start:',start)
    loops = list(nx.simple_cycles(reduced_graph,length_bound=50))
    print('end:',datetime.now())
    print('duration:',datetime.now()-start)
    lengths_loops = np.array([len(l) for l in loops])
    #create a key for the loops and their types
    loop_types = np.array([classify_loop_type(l,new_bonds) for l in loops])
    primary_lengths = lengths_loops[np.where(loop_types == "primary")[0]]
    secondary_lengths = lengths_loops[np.where(loop_types == "secondary")[0]]
    tertiary_lengths = lengths_loops[np.where(loop_types == "tertiary")[0]]
    quaternary_lengths = lengths_loops[np.where(loop_types == "quaternary+")[0]]

    if (testing):
        out_frame = frame
        out_frame.particles.velocity[:] = -1
        for strand in dangling_end_strands:
            for bead in strand:
                out_frame.particles.velocity[bead] = [0,0,0]
        for i,loop in enumerate(loops):
            match loop_types[i]:
                case "primary":
                    type = 1
                case "secondary":
                    type = 2
                case "tertiary":
                    type = 3
                case "quaternary+":
                    type = 4
            for bead in loop:
                out_frame.particles.velocity[bead] = [type,0,0]
        output_gsd = gsd.hoomd.open(jdir + f'/test_{os.path.basename(__file__)}'.replace('.py','.gsd'), 'w')
        output_gsd.append(out_frame)

    unique_primary_lengths,count = np.unique(primary_lengths,return_counts=True)
    primary_data = list(zip(unique_primary_lengths,count))
    print(primary_data)
    unique_secondary_lengths,count = np.unique(secondary_lengths,return_counts=True)
    secondary_data = list(zip(unique_secondary_lengths,count))
    print(secondary_data)
    unique_tertiary_lengths,count = np.unique(tertiary_lengths,return_counts=True)
    tertiary_data = list(zip(unique_tertiary_lengths,count))
    print(tertiary_data)
    unique_quaternary_lengths,count = np.unique(quaternary_lengths,return_counts=True)
    quaternary_data = list(zip(unique_quaternary_lengths,count))
    print(quaternary_data)

    # sizes_of_loops,count = np.unique(lengths_loops,return_counts=True)
    # loop_count = list(zip(sizes_of_loops,count))
    # print(loop_count)
    with open(jdir + "/loop_counts.txt", "w") as f:
        f.write("dangling strand data\n")
        for x in dangling_end_data:
            line = str(x[0]) + " " + str(x[1]) + "\n"
            f.write(line)

        f.write("primary loop_size data\n")
        for x in primary_data:
            line = str(x[0]) + " " + str(x[1]) + "\n"
            f.write(line)

        f.write("secondary loop_size data\n")
        for x in secondary_data:
            line = str(x[0]) + " " + str(x[1]) + "\n"
            f.write(line)

        f.write("tertiary loop_size data\n")
        for x in tertiary_data:
            line = str(x[0]) + " " + str(x[1]) + "\n"
            f.write(line)

        f.write("quaternary loop_size data\n")
        for x in quaternary_data:
            line = str(x[0]) + " " + str(x[1]) + "\n"
            f.write(line)

def xlink_rdf_analysis(job_id):
    print(
    '''
    This function is intended to calculate the rdf of crosslinks in the system.
    However, running this on only 5 frames results in a lot of noise. Thus it needs 
    multiple frames to be run on. For the RDF analysis, there cannot be any new crosslinks
    being formed over the course of this trajectory. This analysis may not be reliable
    enough without additional simulation.
    ''')
    testing = True

    pio.kaleido.scope.mathjax = None

    project = signac.get_project()
    direc = project.fn('') + 'workspace/'

    # for job_id in os.listdir(direc):
    jdir = direc + job_id
    input_file = direc + job_id + "/polymerize.gsd"
    if not os.path.isfile(input_file):
        raise FileNotFoundError("File not found")
        exit()
    # open the polymerized gsd file's trajectory (list of frames)
    trajectory = gsd.hoomd.open(input_file)
    frame = trajectory[-1]
    dummy_id = len(frame.bonds.types)-1
    bonds = frame.bonds.group[frame.bonds.typeid!=dummy_id]

    # make a graph - each bond is an edge
    G = nx.Graph()
    G.add_edges_from(bonds)
    # find all disconnected/connected sub-networks
    # sort the list of networks by length
    Gcc = sorted(nx.connected_components(G), key=len, reverse=True)
        
    # strand lengths (in entire system. If you only want the ones in the gel, you need to look at Gcc[0] only)
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
    #get only the crosslink_beads
    ene_crosslink_beads = [x for  x in G2.nodes() if G2.degree(x) == 3]
    thiol_crosslink_beads = [x for  x in G2.nodes() if G2.degree(x) == 4]
    print("# of ene crosslinks:",len(ene_crosslink_beads))
    print("# of thiol crosslinks:",len(thiol_crosslink_beads))
    crosslinks_only_traj = []
    # make the new frames
    for frame in trajectory[-5:-1]:
        crosslinks_only_frame = frame
        crosslinks_only_frame.particles.position = frame.particles.position[ene_crosslink_beads + thiol_crosslink_beads]
        crosslinks_only_frame.particles.typeid = frame.particles.typeid[ene_crosslink_beads + thiol_crosslink_beads]
        crosslinks_only_frame.particles.orientation = frame.particles.orientation[ene_crosslink_beads + thiol_crosslink_beads]
        crosslinks_only_frame.particles.mass = frame.particles.mass[ene_crosslink_beads + thiol_crosslink_beads]
        crosslinks_only_frame.particles.charge = frame.particles.charge[ene_crosslink_beads + thiol_crosslink_beads]
        crosslinks_only_frame.particles.diameter = frame.particles.diameter[ene_crosslink_beads + thiol_crosslink_beads]
        crosslinks_only_frame.particles.body = frame.particles.body[ene_crosslink_beads + thiol_crosslink_beads]
        crosslinks_only_frame.particles.moment_inertia = frame.particles.moment_inertia[ene_crosslink_beads + thiol_crosslink_beads]
        crosslinks_only_frame.particles.velocity = frame.particles.velocity[ene_crosslink_beads + thiol_crosslink_beads]
        crosslinks_only_frame.particles.angmom = frame.particles.angmom[ene_crosslink_beads + thiol_crosslink_beads]
        crosslinks_only_frame.particles.image = frame.particles.image[ene_crosslink_beads + thiol_crosslink_beads]
        crosslinks_only_frame.particles.N = len(ene_crosslink_beads) + len(thiol_crosslink_beads)

        crosslinks_only_traj.append(crosslinks_only_frame)
    #get the rdfs
    # 'C' = tetrafunctional thiol crosslinker center, 'Carbon' = chain growthed ene
    rdf_ene, norm_ene = intermolecular_rdf(crosslinks_only_traj, 
                                           A_name='Carbon', B_name='Carbon', 
                                           r_max=2.0, r_min=0.25, bins=100, 
                                           exclude_bonded=False,)
    rdf_thiol, norm_thiol = intermolecular_rdf(crosslinks_only_traj, 
                                               A_name='C', B_name='C', 
                                               r_max=2.0, r_min=0.25, bins=100, 
                                               exclude_bonded=False,)
    rdf_thiol_ene, norm_thiol_ene = intermolecular_rdf(crosslinks_only_traj, 
                                                       A_name='C', B_name='Carbon', 
                                                       r_max=2.0, r_min=0.25, bins=100, 
                                                       exclude_bonded=False,) 

    # save the rdf data into txts
    with open(jdir + "/ene_ene_rdf.txt", "w") as f:
        f.write("r g(r)\n")
        for r, g in zip(rdf_ene.bin_centers, rdf_ene.rdf*norm_ene):
            f.write(f"{r} {g}\n")
    with open(jdir + "/thiol_thiol_rdf.txt", "w") as f:
        f.write("r g(r)\n")
        for r, g in zip(rdf_thiol.bin_centers, rdf_thiol.rdf*norm_thiol):
            f.write(f"{r} {g}\n")
    with open(jdir + "/thiol_ene_rdf.txt", "w") as f:
        f.write("r g(r)\n")
        for r, g in zip(rdf_thiol_ene.bin_centers, rdf_thiol_ene.rdf*norm_thiol_ene):
            f.write(f"{r} {g}\n")

    if testing:
        # make the plots
        fig = make_subplots(rows=1, cols=3, subplot_titles=("Crosslink Ene-Ene RDF", "Crosslink Thiol-Thiol RDF", "Crosslink Thiol-Ene RDF"))
        fig.add_trace(go.Scatter(x=rdf_ene.bin_centers, y=rdf_ene.rdf*norm_ene, mode='lines', name='Ene-Ene RDF'), row=1, col=1)
        fig.add_trace(go.Scatter(x=rdf_thiol.bin_centers, y=rdf_thiol.rdf*norm_thiol, mode='lines', name='Thiol-Thiol RDF'), row=1, col=2)
        fig.add_trace(go.Scatter(x=rdf_thiol_ene.bin_centers, y=rdf_thiol_ene.rdf*norm_thiol_ene, mode='lines', name='Thiol-Ene RDF'), row=1, col=3)
        fig.show()
        fig.write_image("./test_network_properties/rdfs.jpg")


def main(job_id):
    xlink_rdf_analysis(job_id)
    exit()
    defect_analysis(job_id)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 network_properties.py job_id")
        sys.exit(1)
    main(sys.argv[1])