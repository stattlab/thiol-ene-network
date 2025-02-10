#!/opt/miniconda3/bin/python
"""Define the project's workflow logic and operation functions.

Execute this script directly from the command line, to view your project's
status, execute operations and submit them to a cluster. See also:

    $ python project.py --help
"""
import flow
from flow import FlowProject
import os
import subprocess
import numpy as np

class MyProject(FlowProject):
    pass

class CampusCluster(flow.environment.DefaultSlurmEnvironment):
    hostname_pattern = r".*.campuscluster\.illinois\.edu$"
    #template = "campuscluster.sh"

# Labels
@MyProject.label
def equilibrated(job):
    return job.isfile("equi.gsd")

@MyProject.label
def reacted(job):
    # return job.doc['reacted_monomers']>0.9
    try:
        conversion_data = np.genfromtxt(job.fn('trajectory_conversion.txt'), comments="#", delimiter=" ")
    except FileNotFoundError:
        return False
    cutoff = 0.00
    if float(conversion_data[-1,1]) - float(conversion_data[-50,1]) < cutoff:
        return True
    else:
        return False

@MyProject.label
def reacting(job):
    return job.isfile("polymerize.gsd")

@MyProject.label
def deformed(job):
    return job.isfile("deform.gsd") and os.path.getsize(job.fn('deform.log')) > 100

#-----------------------
# Analysis labels
#-----------------------

# def percolation_analyzed(job):
#     return job.isfile('percolation_data.txt')


#-----------------------
# Simulation operations
#-----------------------
@MyProject.post(equilibrated)
@MyProject.operation
def equilibrate(job):
    from scripts.simulate import Simulator
    # if job.sp["polymerization_method"] == "cpu_local_snapshot":
    #     from scripts.simulate import Simulator
    # elif job.sp["polymerization_method"] == "custom_action_GPU" or job.sp["polymerization_method"] == "custom_action_CPU" or job.sp["polymerization_method"] == "custom_action_GPU_bulk":
    #     from scripts.customAction.simulate import Simulator
    # else:
    #     print("ERROR polymerization_method not recognized")
    #     exit(2)
    sinit = Simulator(job)
    sinit.equilibrate()
    print("equilibrated",job.id)


@MyProject.pre(equilibrated)
@MyProject.post(reacted)
@MyProject.operation
def polymerize(job):
    from scripts.simulate import Simulator
    # if job.sp["polymerization_method"] == "cpu_local_snapshot":
    #     from scripts.simulate import Simulator
    # elif job.sp["polymerization_method"] == "custom_action_GPU" or job.sp["polymerization_method"] == "custom_action_CPU" or job.sp["polymerization_method"] == "custom_action_GPU_bulk":
    #     from scripts.customAction.simulate import Simulator
    # else:
    #     print("ERROR polymerization_method not recognized")
    #     exit(2)
    sinit = Simulator(job)
    sinit.polymerize()
    print("reacted",job.id)
    print("now conducting trajectory conversion analysis")
    os.system('python3 ./scripts/trajectory_conversion.py')
    print('analyzed the trajectory conversion of ', job.id)

@MyProject.pre(reacted)
@MyProject.post(deformed)
@MyProject.operation
def deform(job):
    import scripts.deform
    scripts.deform.main(job.fn('polymerize.gsd'),-1)
    print('deformed: ',job.id)

#-----------------------
# Analysis operations
#-----------------------
# @MyProject.pre(reacting)
# @MyProject.post(trajectory_conversion_analyzed)
# @MyProject.operation
# def analyze_trajectory_conversion(job):
#     os.system('python3 ./scripts/trajectory_conversion.py')





if __name__ == "__main__":
    MyProject().main()
