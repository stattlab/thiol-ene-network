import sys, os, re, json, io, itertools
import numpy as np
import signac
import gsd, gsd.hoomd 
from collections import defaultdict
from collections import OrderedDict
from collections import Counter


"""
frame.bonds.group gives a np.array that is of form [[a,b],[c,d],[e,f],[a,g],...] 
where a is bonded to b, c is bonded to d, etc.

atom_bonds finds all atoms to which a given atom 'a' is bonded to
args: frame.bonds.group

returns: dictionary of atoms and the atoms they are bonded to
{atom_no: [bonded1, bonded2, ...], ...}
"""
def atom_bonds(bonds):
    # create a default dictionary with an element returning a list if no neighbors
    neighbors = defaultdict(list) 
    # add every instance of an atom in a bond to the dictionary
    for each in bonds:
        for item in each:
            neighbors[item].extend(each)
            # remove the atom itself from the list of atoms bonded to it
            neighbors[item].remove(item)
    return neighbors


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
'''
# make a default figure box 
def make_plotly_fig(xaxis="", xaxis_range=None, yaxis="", yaxis_range=None):
    fig = go.Figure()
    fig.update_layout(
        font_family="Avenir Medium",
        font_color="black",
    )
    fig.update_layout(
        plot_bgcolor='white'
    )
    fig.update_layout(showlegend=False,
                      width=500, 
                      height =500,
                      xaxis_title=xaxis,
                      yaxis_title=yaxis)
    
    fig.update_xaxes(linecolor='black', mirror=True,
                     zeroline=False,
                     title_font_size=20,
                     tickfont_size = 18,
                     ticks='inside'
                    )
    
    fig.update_yaxes(mirror=True,ticks='inside',linecolor='black', zeroline=False,
                     title_font_size=20,
                     tickfont_size = 20,
                    )
    if xaxis_range != None:
        fig.update_xaxes(range=xaxis_range)
    if yaxis_range != None:
        fig.update_yaxes(range=yaxis_range)
    return fig
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream

def analyze_youngs_modulus(jobid):
    project = signac.get_project()
    job = project.open_job(id=jobid)
    data = np.genfromtxt(job.fn('stress_strain.txt'))
    strain = data[:,0]
    stress = data[:,1]
    #---------------------------------------------------------------------
    #                    Calculate and save Young's Modulus
    #---------------------------------------------------------------------
    # find the linear region of the stress-strain curve up to 2% (https://www.sciencedirect.com/science/article/pii/S0032386112007318)
    linear_strains = strain[:np.argmax(strain > 0.02)]
    linear_stresses = stress[:np.argmax(strain > 0.02)]
    print(f"Calculating Young's Modulus with {len(linear_strains)} points")
    # find the slope of the linear region
    slope, intercept = np.polyfit(linear_strains, linear_stresses, 1)
    # save the Young's modulus
    with open(job.fn('youngs_modulus.txt'), "w") as f:
        f.write(str(slope))

def analyze_deform(jobid):
    project = signac.get_project()
    job = project.open_job(id=jobid)
    # get and trajectory and make sure the first frame checks out
    trajectory = gsd.hoomd.open(job.fn('pressures.gsd'))
    logs = gsd.hoomd.read_log(job.fn('pressures.gsd'))
    pressure_tensors = logs["log/md/compute/ThermodynamicQuantities/pressure_tensor"]
    data = np.genfromtxt(job.fn('deform.log'))
    Lx_arr = data[1:,7]
    Lx0 = Lx_arr[0]

    strain = []
    true_stress_deviatoric = []
    true_stress_straight = []
    # iterate through each frame
    for i,frame in enumerate(trajectory):
        strain.append(np.log(Lx_arr[i]/Lx0))#calculate true strain
        #two ways of calculating true stress
        #1
        true_stress_straight.append(-1*pressure_tensors[i][0])
        #2
        hydrostaticPressure = np.add(pressure_tensors[i][0],np.add(pressure_tensors[i][3],pressure_tensors[i][5]))/3
        deviatoricPressure = np.subtract(pressure_tensors[i][0],hydrostaticPressure)
        true_stress_deviatoric.append(deviatoricPressure)
    unique_strains, indices = np.unique(strain, return_index=True)
    unique_true_stress_straight = []
    for i in range(len(indices)-1):
        unique_true_stress_straight.append(np.average(true_stress_straight[indices[i]:indices[i+1]]))
    unique_true_stress_straight.append(np.average(true_stress_straight[indices[-1]:]))

    # window averaging
    N = 5
    avg_unique_strains = np.convolve(unique_strains,np.ones(N)/N,mode='valid')
    avg_unique_true_stress_straight = np.convolve(unique_true_stress_straight,np.ones(N)/N,mode='valid')

    # Save the stress-strain curve
    with open(job.fn('stress_strain.txt'), "w") as f:
        for i in range(len(avg_unique_strains)):
            f.write(f"{avg_unique_strains[i]} {avg_unique_true_stress_straight[i]}\n")


def main():
    plt.rcParams["font.family"] = "Avenir"
    color = iter(cm.rainbow(np.linspace(0, 1, 11)))
    fig, ax = plt.subplots(1,1,sharey=False)

    '''
    WHAT TO EDIT**:

    ------------------------------------------------------

    '''
=======
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
'''

def analyze_modulus(jobid):
    project = signac.get_project()
    job = project.open_job(id=jobid)
    data = np.genfromtxt(job.fn('stress_strain.txt'))
    strain = data[:,0]
    stress = data[:,1]
    #---------------------------------------------------------------------
    #                    Calculate and save Young's Modulus
    #---------------------------------------------------------------------
    # find the linear region of the stress-strain curve up to 2% (https://www.sciencedirect.com/science/article/pii/S0032386112007318)
    linear_strains = strain[:np.argmax(strain > 0.02)]
    linear_stresses = stress[:np.argmax(strain > 0.02)]
    print(f"Calculating Young's Modulus with {len(linear_strains)} points")
    # find the slope of the linear region
    slope, intercept = np.polyfit(linear_strains, linear_stresses, 1)
    # save the Young's modulus
    with open(job.fn('youngs_modulus.txt'), "w") as f:
        f.write(str(slope))

def analyze_stress_strain(jobid):
    project = signac.get_project()
    job = project.open_job(id=jobid)
    # get and trajectory and make sure the first frame checks out
    trajectory = gsd.hoomd.open(job.fn('pressures.gsd'))
    logs = gsd.hoomd.read_log(job.fn('pressures.gsd'))
    pressure_tensors = logs["log/md/compute/ThermodynamicQuantities/pressure_tensor"]
    data = np.genfromtxt(job.fn('deform.log'))
    Lx_arr = data[1:,7]
    Lx0 = Lx_arr[0]

    strain = []
    true_stress_deviatoric = []
    true_stress_straight = []
    # iterate through each frame
    for i,frame in enumerate(trajectory):
        strain.append(np.log(Lx_arr[i]/Lx0))#calculate true strain
        #two ways of calculating true stress
        #1
        true_stress_straight.append(-1*pressure_tensors[i][0])
        #2
        hydrostaticPressure = np.add(pressure_tensors[i][0],np.add(pressure_tensors[i][3],pressure_tensors[i][5]))/3
        deviatoricPressure = np.subtract(pressure_tensors[i][0],hydrostaticPressure)
        true_stress_deviatoric.append(deviatoricPressure)
    unique_strains, indices = np.unique(strain, return_index=True)
    unique_true_stress_straight = []
    for i in range(len(indices)-1):
        unique_true_stress_straight.append(np.average(true_stress_straight[indices[i]:indices[i+1]]))
    unique_true_stress_straight.append(np.average(true_stress_straight[indices[-1]:]))

    # window averaging
    N = 5
    avg_unique_strains = np.convolve(unique_strains,np.ones(N)/N,mode='valid')
    avg_unique_true_stress_straight = np.convolve(unique_true_stress_straight,np.ones(N)/N,mode='valid')

    # Save the stress-strain curve
    with open(job.fn('stress_strain.txt'), "w") as f:
        for i in range(len(avg_unique_strains)):
            f.write(f"{avg_unique_strains[i]} {avg_unique_true_stress_straight[i]}\n")

'''
def main():
    plt.rcParams["font.family"] = "Avenir"
    color = iter(cm.rainbow(np.linspace(0, 1, 11)))
    fig, ax = plt.subplots(1,1,sharey=False)
<<<<<<< Updated upstream
<<<<<<< Updated upstream
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
    # our signac project
    project = signac.get_project()
    # look through each POLYMERIZED job in the project

    # group_statepoint = {"replica_index": 8, "density": 0.9, "temperature": 0.9, "crosslinker_percent": 50, 
    #                     "N_monomers": 500, "monomer_size": 0, "extender_size": 0, "radical_percent": 1.0, 
    #                     "chain_side_reaction_probability": 0, "chain_transfer_probability": 0.5, 
    #                     "thiol_reaction_probability": 0.5, "polymerize_period": 100, "r_cut_reaction": 1.1, 
    #                     "angle_constant": 5.0}
    # # which group to plot
    # group_plot_name = "0_mon_500_cg"
    # group_plot_dir = "./plots/" + group_plot_name + "/"
    # if not os.path.isdir(group_plot_dir):
    #     os.mkdir(group_plot_dir)
    # cg_only = True


<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
    '''
    END**
    ------------------------------------------------------
    '''

    for job in project: 
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
        # if the job has been deformed
        if job.isfile('pressures.gsd'):
            with open(job.fn('signac_statepoint.json')) as f:
                statepoint = json.load(f)

            # try:
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
            ''' initial look '''
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
            # get and trajectory and make sure the first frame checks out
            trajectory = gsd.hoomd.open(job.fn('pressures.gsd'))

            logs = gsd.hoomd.read_log(job.fn('pressures.gsd'))
            pressure_tensors = logs["log/md/compute/ThermodynamicQuantities/pressure_tensor"]

            data = np.genfromtxt(job.fn('deform.log'))
            # print(data)
            Lx_arr = data[1:,7]
            print(len(Lx_arr))
            print(len(trajectory))
            # print(np.shape(Lx_arr))

            #analyze
            strain = []
            true_stress_deviatoric = []
            true_stress_straight = []

            # init_box = trajectory[0].configuration.box
            Lx0 = Lx_arr[0]

            # iterate through each frame
            for i,frame in enumerate(trajectory):
                strain.append(np.log(Lx_arr[i]/Lx0))#calculate true strain
                
                #two ways of calculating true stress
                #1
                true_stress_straight.append(-1*pressure_tensors[i][0])
                #2
                hydrostaticPressure = np.add(pressure_tensors[i][0],np.add(pressure_tensors[i][3],pressure_tensors[i][5]))/3
                deviatoricPressure = np.subtract(pressure_tensors[i][0],hydrostaticPressure)
                true_stress_deviatoric.append(deviatoricPressure)

            unique_strains, indices = np.unique(strain, return_index=True)
            unique_true_stress_straight = []
            for i in range(len(indices)-1):
                unique_true_stress_straight.append(np.average(true_stress_straight[indices[i]:indices[i+1]]))

            unique_true_stress_straight.append(np.average(true_stress_straight[indices[-1]:]))
            
            ### plots ------------------------------------------
            plot_dir = "./workspace/" + str(job) + "/plots/"
            if not os.path.isdir(plot_dir):
                os.mkdir(plot_dir)

            # window averaging
            N = 5
            avg_unique_strains = np.convolve(unique_strains,np.ones(N)/N,mode='valid')
            avg_unique_true_stress_straight = np.convolve(unique_true_stress_straight,np.ones(N)/N,mode='valid')
            ax.plot(avg_unique_strains,
                    avg_unique_true_stress_straight,
                    linewidth=2,label=rf'{job.id}')#$\sigma_{xx}$')
            plt.legend()
            # ax.plot(unique_strains,unique_true_stress_straight)
            # ax[1].plot(strain,true_stress_deviatoric)
            # plt.savefig(plot_dir + 'stress_strain.png')

            plt.savefig(plot_dir + 'stress_strain.jpg')


            #---------------------------------------------------------------------
            #                    Calculate and save Young's Modulus
            #---------------------------------------------------------------------
            # find the linear region of the stress-strain curve up to 2% (https://www.sciencedirect.com/science/article/pii/S0032386112007318)
            linear_strains = avg_unique_strains[:np.argmax(avg_unique_strains > 0.02)]
            linear_stresses = avg_unique_true_stress_straight[:np.argmax(avg_unique_strains > 0.02)]
            print(f"Calculating Young's Modulus with {len(linear_strains)} points")
            # find the slope of the linear region
            slope, intercept = np.polyfit(linear_strains, linear_stresses, 1)
            # save the Young's modulus
            with open(job.fn('youngs_modulus.txt'), "w") as f:
                f.write(str(slope))
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream

if __name__ == '__main__':
    main()
=======
'''
if __name__ == '__main__':
    main()
>>>>>>> Stashed changes
=======
'''
if __name__ == '__main__':
    main()
>>>>>>> Stashed changes
=======
'''
if __name__ == '__main__':
    main()
>>>>>>> Stashed changes
