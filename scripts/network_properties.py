import sys, os, re, json, io, itertools
import numpy as np
import signac
import gsd, gsd.hoomd 
from collections import defaultdict
from collections import OrderedDict
from collections import Counter
from datetime import datetime

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

def classify_loop_type_simple(cycle):
    #cycle is the list of nodes in the loop
    #for this classification to be correct, these nodes need to only include crosslinkers
    match len(cycle):
        case 1:
            return "primary"
        case 2:
            return "secondary"
        case 3:
            return "tertiary"
        case 4:
            return "quaternary+"
        case _:
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

def reduce_network_dangling_ends(in_graph):
    '''
    Does one iteration of removing dangling ends from a network.
    Dangling ends are defined as nodes with only one edge, and must be evaluated recursively
        to account for branching dangling ends. (see remove_dangling_ends())
    Args:
        in_graph: networkx.Graph instance
    Returns:
        out_graph: networkx.Graph instance with one iteration of dangling ends removed
        dangling_ends: list of beads that were a part of dangling ends removed from the network
    '''
    dangling_ends = []
    out_graph = in_graph.copy()
    for i in in_graph.nodes():
        if in_graph.degree(i) == 1:
            out_graph.remove_node(i)
            dangling_ends.append(i)
    return out_graph, dangling_ends

#returns a MULTIGRAPH
def remove_extenders(in_graph):
    ''''
    Recursively removes all nodes with degree 2, replacing them with a single edge 
        connecting the two neighbors. Essentially generates a CROSSLINK-ONLY graph.
    Args:
        in_graph: networkx.Graph instance
    Returns:
        out_graph: networkx.MultiGraph instance with only crosslinks
    '''
    # out_graph = in_graph.copy().to_multigraph()
    out_graph = nx.MultiGraph()
    out_graph.add_edges_from(in_graph.edges(),weight=0)
    out_graph.add_nodes_from(in_graph.nodes(),weight=0)
    strands,crosslinkers = cut_into_strands(in_graph)
    for strand in strands:
        #get the nodes that are connected to crosslinkers
        #get the crosslinkers that the end nodes are connected to
        end_nodes = []
        strand_crosslinkers = []
        for node in strand:
            neighbors = list(in_graph.neighbors(node))
            for neighbor in neighbors:
                if neighbor not in strand:
                    strand_crosslinkers.append(neighbor)
            if len(neighbors) == 1:
                end_nodes.append(node)
        # remove all the nodes in the strand, replacing them with an edge connecting the crosslinkers
        out_graph.remove_nodes_from(strand)
        if len(strand_crosslinkers) == 2:
            # print("adding edge between:",strand_crosslinkers[0],strand_crosslinkers[1])
            out_graph.add_edge(strand_crosslinkers[0],strand_crosslinkers[1],weight=len(strand))
        elif len(strand_crosslinkers) == 1:
            # print("adding edge between:",strand_crosslinkers[0],end_nodes[0])
            out_graph.add_edge(strand_crosslinkers[0],end_nodes[0],weight=len(strand))
        else:
            print(f"Strand {strand} has {len(strand_crosslinkers)} ends/crosslinkers: {strand_crosslinkers}")
            raise ValueError("Strand has more than 2 ends/crosslinkers")
    return out_graph

def cut_into_strands(in_graph):
    '''
    Takes a graph and cuts it into strands which end right before any crosslinkers with 
        functionality 3 or higher

    Args:
        in_graph: networkx.Graph instance
    Returns:
        Gcc_new: list of networkx.Graph instances, each representing a strand sorted largest to smallest
        crosslink_beads: list of beads that were removed from the network
    '''
    Gnew = in_graph.copy()
    crosslink_beads = [x for  x in in_graph.nodes() if in_graph.degree(x) >= 3]
    
    for x in crosslink_beads:
        Gnew.remove_node(x)

    # strands
    Gcc_new = sorted(nx.connected_components(Gnew), key=len, reverse=True)
    return Gcc_new,crosslink_beads

def remove_dangling_ends(in_graph):
    '''
    Recursively removes dangling ends from a network. 
    Dangling ends are defined as nodes with only one edge, and must be evaluated recursively
        to account for branching dangling ends.
    Args:
        in_graph: unmodified networkx.Graph instance
    Returns:
        new_graph: networkx.Graph instance with dangling ends removed
        dangling_ends: list of beads that were a part of dangling ends removed from the network
    '''
    dangling_ends = []
    new_graph = nx.Graph()
    last_graph = in_graph.copy()
    for i in range(10000):
        new_graph, new_dangling_ends = reduce_network_dangling_ends(last_graph)
        dangling_ends.extend(new_dangling_ends)
        if new_graph.edges() == last_graph.edges():
            break
        else:
            last_graph = new_graph.copy()
    else:
        raise ValueError("Dangling ends not removed after 10,000 iterations. Consider increasing the limit, or check for infinite recursion.")
    return new_graph, dangling_ends

# To contract chains of neighbouring vertices with degree 2 into one hypernode:
#https://stackoverflow.com/questions/52313551/graph-reduction

def is_in_2D(list, element):
    """
    Checks if an element exists in a two-dimensional list.

    Args:
        list: The two-dimensional list to search within.
        element: The element to search for.

    Returns:
        The index of the sublist containing the element, or -1 if not found.
    """
    sublist_index = 0
    for sublist in list:
        if element in sublist:
            return sublist_index
        sublist_index += 1
    return -1 #Not found

def categorize_loop_edges(loops,crosslink_graph):
    '''
    Categorizes the edges in each loop into primary, secondary, tertiary, and quaternary+ loops.
    Assumes that the crosslink_graph is a graph of only crosslinks, where the nodes in 
        between have been replaced with a single edge.
    Primary loops are loops that contain only one edge.
    Secondary loops are loops that contain two edges.
    Tertiary loops are loops that contain three edges.
    Quaternary+ loops are loops that contain four or more edges.

    If an edge is shared between two types of loop, it is categorized as the lower degree type.
    
    Args:
        loops: list of lists of nodes in each loop
        crosslink_graph: networkx.Graph instance of ONLY CROSSLINKS (see function remove_extenders)
    Returns:
        primary_edges: list of lists of edges in primary loops
        secondary_edges: list of lists of edges in secondary loops
        tertiary_edges: list of lists of edges in tertiary loops
        quaternary_edges: list of lists of edges in quaternary+ loops
    '''
    #graph needs to be the graph of only crosslinks
    loop_edges = np.empty((len(loops),),dtype=object)
    loop_types = np.array([classify_loop_type_simple(l) for l in loops])
    #returns lists of the edges in each loop type
    primary_edges = []
    secondary_edges = []
    tertiary_edges = []
    quaternary_edges = []

    #Compile lists of the edges in each loop
    for i,loop in enumerate(loops):
        loop_edges[i] = []
        if len(loop) == 1:
            loop_edges[i].append(crosslink_graph.get_edge_data(u=loop[0],v=loop[0]))
        else:
            for j in range(len(loop)):
                for k in range(j+1,len(loop)):
                    bead0=loop[j]
                    bead1=loop[k]
                    loop_edges[i].append(crosslink_graph.get_edge_data(u=bead0,v=bead1))
        if loop_types[i] == "primary":
            primary_edges.append(loop_edges[i])
        elif loop_types[i] == "secondary":
            secondary_edges.append(loop_edges[i])
        elif loop_types[i] == "tertiary":
            tertiary_edges.append(loop_edges[i])
        elif loop_types[i] == "quaternary+":
            quaternary_edges.append(loop_edges[i])

    #Remove double counting of edges
    for i in range(len(primary_edges)):
        test_index = is_in_2D(element=primary_edges[i], list=secondary_edges)
        if (test_index != -1):
            secondary_edges[test_index].remove(primary_edges[i])
        test_index = is_in_2D(element=primary_edges[i], list=tertiary_edges)
        if (test_index != -1):
            tertiary_edges[test_index].remove(primary_edges[i])
        test_index = is_in_2D(element=primary_edges[i], list=quaternary_edges)
        if (test_index != -1):
            quaternary_edges[test_index].remove(primary_edges[i])
    for i in range(len(secondary_edges)):
        for j in range(len(secondary_edges[i])):
            test_index = is_in_2D(element=secondary_edges[i][j], list=tertiary_edges)
            if (test_index != -1):
                tertiary_edges[test_index].remove(secondary_edges[i][j])
            test_index = is_in_2D(element=secondary_edges[i][j], list=quaternary_edges)
            if (test_index != -1):
                quaternary_edges[test_index].remove(secondary_edges[i][j])
        # if (is_in_2D(element=secondary_edges[i], list=tertiary_edges):
        #     tertiary_edges.remove(secondary_edges[i])
        # if (is_in_2D(element=secondary_edges[i], list=quaternary_edges):
        #     quaternary_edges.remove(secondary_edges[i])
    for i in range(len(tertiary_edges)):
        test_index = is_in_2D(element=tertiary_edges[i], list=quaternary_edges)
        if test_index != -1:
            quaternary_edges[test_index].remove(tertiary_edges[i])

    for i in range(len(quaternary_edges)):
        while None in quaternary_edges[i]:
            quaternary_edges[i].remove(None)
        # if quaternary_edges[i] == []:
        #     quaternary_edges.pop(i)
        # if is_in_2D(element=None,list=quaternary_edges) != -1:
        #     print(quaternary_edges[index])
        #     quaternary_edges.pop(index)

    return primary_edges, secondary_edges, tertiary_edges, quaternary_edges

def get_length_categorized_loops(primary_edges,secondary_edges,tertiary_edges,quaternary_edges):
    '''
    Returns the length of each loop in each category of loop type.
    Args:
        primary_edges: list of lists of edges in primary loops
        secondary_edges: list of lists of edges in secondary loops
        tertiary_edges: list of lists of edges in tertiary loops
        quaternary_edges: list of lists of edges in quaternary+ loops
    Returns:
        primary_lengths: list of lengths of primary loops
        secondary_lengths: list of lengths of secondary loops
        tertiary_lengths: list of lengths of tertiary loops
        quaternary_lengths: list of lengths of quaternary+ loops
    '''
    #graph needs to be the graph of only crosslinks
    primary_lengths = np.zeros(len(primary_edges))
    secondary_lengths = np.zeros(len(secondary_edges))
    tertiary_lengths = np.zeros(len(tertiary_edges))
    quaternary_lengths = np.zeros(len(quaternary_edges))

    for i,loop_edges in enumerate(primary_edges):
        # print("primary edges:",edges)
        tmp = 0
        for edge in loop_edges:
            for k in range(len(edge)):
                edge_weight = edge[k]['weight']
                tmp += edge_weight
        primary_lengths[i] = tmp

    for i,loop_edges in enumerate(secondary_edges):
        # print("secondary edges:",edges)
        tmp = 0
        for edge in loop_edges:
            for k in range(len(edge)):
                edge_weight = edge[k]['weight']
                tmp += edge_weight
        secondary_lengths[i] = tmp

    for i,loop_edges in enumerate(tertiary_edges):
        # print("tertiary edges:",edges)
        tmp = 0
        for edge in loop_edges:
            for k in range(len(edge)):
                edge_weight = edge[k]['weight']
                tmp += edge_weight
        tertiary_lengths[i] = tmp

    for i,loop_edges in enumerate(quaternary_edges):
        tmp = 0
        for edge in loop_edges:
            for k in range(len(edge)):
                edge_weight = edge[k]['weight']
                tmp += edge_weight
        quaternary_lengths[i] = tmp

    return primary_lengths, secondary_lengths, tertiary_lengths, quaternary_lengths

#---------------------------------------------------------------------------------------
# Analysis Functions
#---------------------------------------------------------------------------------------

#Dangling end lengths do not include the beads within the gel
#
def defect_analysis(job_id,testing=False):
    
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


    Ntotal = frame.particles.N

    if testing:
        jdir = "./test_network_properties"
        bonds = [[0,1],
                 [1,1.5],[1.5,2],[1,1.6],[1.6,2],[11,2.11],[2,2.11],
                 [2,3],[3,4],[4,1],[0,5],[5,6],[6,7],[5,8],[8,9],
                    [2,10],[2,11],[10,12],[12,11],[12,13],[12,14],[13,15],[15,14],[0,16],[0,11],[6,17],[1,1]]
        print("bonds:",bonds)
        G = nx.Graph()
        G.add_edges_from(bonds)
        G2 = nx.Graph()
        G2.add_edges_from(bonds)
        frame = gsd.hoomd.Frame()
        frame.particles.position = np.array([[0,0,0],[1,0,0],[1.5,0,0],[1.6,0,0],[2,0,0],[2.11,0,0],[2,1,0],[1,1,0],[-1,0,0],[-2,0,0],[-3,0,0],[-1,1,0],[-1,2,0],
                                             [3,0,0],[2,1,0],[3,1,0],[4,1,0],[3,2,0],[4,2,0],[-1,0,0],[-2,-1,0]])
        frame.particles.velocity = np.zeros((len(frame.particles.position),3))
        frame.particles.N = len(frame.particles.position)
        frame.bonds.N = len(bonds)
        frame.bonds.group = np.array(bonds)
        frame.bonds.typeid = np.zeros(len(bonds))
        frame.bonds.types = ['bond']
        new_bonds = bonds

        # print('G2 Nodes:',G2.nodes())
        # print('G2 Edges:',G2.edges())
        # G2 = remove_extenders(G2)
        # print('G2 Nodes:',G2.nodes())
        # print('G2 Edges:',G2.edges())

        '''
        trajectory = gsd.hoomd.open('/home/bj21/simulations/thiol-ene/network_polymerization/workspace/b00fe99f07ce7c0bc08a9cf517888351/polymerize.gsd')

        frame = trajectory[-1]
        dummy_id = len(frame.bonds.types)-1
        bonds = frame.bonds.group[frame.bonds.typeid!=dummy_id]
        '''

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
    else:
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
    
    print("calculating dangling ends")
    dangling_end_strands = []
    dangling_end_lengths = []

    no_danglingEnds_graph, dangling_ends = remove_dangling_ends(G2)
    # print('no_danglingEnds_graph edges', no_danglingEnds_graph.edges())
    # print('no_danglingEnds_graph nodes', no_danglingEnds_graph.nodes())
    dangling_end_graph = G2.copy()
    for x in G2.nodes():
        if x not in dangling_ends:
            dangling_end_graph.remove_node(x)
    # print("dangling ends:",dangling_end_graph.nodes())
    # print("dangling end edges:",dangling_end_graph.edges())

    for dangling_strand in nx.connected_components(dangling_end_graph):
        dangling_strand = list(dangling_strand)
        # print(dangling_strand)
        dangling_end_strands.append(dangling_strand)
        dangling_end_lengths.append(len(dangling_strand))
    unique_dangling_end_lengths,count = np.unique(dangling_end_lengths,return_counts=True)
    dangling_end_data = list(zip(unique_dangling_end_lengths,count))
    # print(dangling_end_data)

    print("calculating loops",  flush=True)
    start = datetime.now()
    crosslink_graph = remove_extenders(no_danglingEnds_graph)
    # print('crosslink_graph edges', crosslink_graph.edges())
    # print('crosslink_graph nodes', crosslink_graph.nodes())
    print('start loop analysis:',start,  flush=True)
    loops = list(nx.simple_cycles(crosslink_graph,length_bound=4),  flush=True)
    print('end loop analysis:',datetime.now(),  flush=True)
    print('duration:',datetime.now()-start,  flush=True)
    # lengths_loops = np.array([len(l) for l in loops])
    loop_types = np.array([classify_loop_type_simple(l) for l in loops])
    primary_edges, secondary_edges, tertiary_edges, quaternary_edges = categorize_loop_edges(loops, crosslink_graph)
    primary_lengths, secondary_lengths, tertiary_lengths, quaternary_lengths = get_length_categorized_loops(primary_edges,secondary_edges,tertiary_edges,quaternary_edges)
    # print('primary edges:',primary_edges)
    # print('secondary edges:',secondary_edges)
    # print('tertiary edges:',tertiary_edges)
    # print('quaternary edges:',quaternary_edges)
    #if testing, print out the ids of the particles participating in loops

    if (testing):
        for i in range(len(loops)):
            print("identified loop:",loops[i]," as type:",loop_types[i])

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

    #Format the loop data as discrete histograms
    unique_primary_lengths,count = np.unique(primary_lengths,return_counts=True)
    primary_data = list(zip(unique_primary_lengths,count))
    print('primary lengths and counts:',primary_data)
    unique_secondary_lengths,count = np.unique(secondary_lengths,return_counts=True)
    secondary_data = list(zip(unique_secondary_lengths,count))
    print('secondary lengths and counts:',secondary_data)
    unique_tertiary_lengths,count = np.unique(tertiary_lengths,return_counts=True)
    tertiary_data = list(zip(unique_tertiary_lengths,count))
    print('tertiary lengths and counts:',tertiary_data)
    unique_quaternary_lengths,count = np.unique(quaternary_lengths,return_counts=True)
    quaternary_data = list(zip(unique_quaternary_lengths,count))
    print('quarternary lengths and counts:', quaternary_data)

    #write the data to a file
    # Here, each line is a tuple of the form (loop_size, count), with header lines 
    # indicating the type of loop data to follow until the next header
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
    G = nx.MultiGraph()
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


def strand_lengths_analysis(job_id):
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
    # calculate gel fraction
    sizes = [len(n) for n in Gcc ]
    total = sum(sizes)
    sizes = sorted(Counter(sizes).items(), key=lambda item: item[0], reverse=True)
    print(sizes[0][0])
    gf = sizes[0][0]/total

    with open(jdir + "/gel_fraction.txt", "w") as f:
        f.write(str(gf))

    new_bonds = []
    for each in bonds:
        for i in each:
            if i in Gcc[0]:
                new_bonds.append(each)
                continue
    # make cluster of biggest graph
    G2 = nx.Graph()
    G2.add_edges_from(new_bonds)
    Gnew = G2.copy()
    # remove everyone that has 3 or more bonds on it, only leaving linear strands
    crosslink_beads = [x for  x in G2.nodes() if G2.degree(x) >= 3]
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
    
    crosslink_density = len(crosslink_beads)/sum([i[0]*i[1] for i in strand_count])
    with open(jdir + "/crosslink_density_1.txt", "w") as f:
        f.write(str(crosslink_density))

    
    




def main(job_id):
    xlink_rdf_analysis(job_id)
    exit()
    defect_analysis(job_id)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 network_properties.py job_id")
        sys.exit(1)
    main(sys.argv[1])
