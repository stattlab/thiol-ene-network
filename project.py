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
    return job.doc['reacted_monomers']>0.9 

@MyProject.label
def reacting(job):
    return job.isfile("polymerize.gsd")


@MyProject.post(equilibrated)
@MyProject.operation
def equilibrate(job):
    from scripts.simulate import Simulator
    sinit = Simulator(job)
    sinit.equilibrate()
    print("equilibrated",job.id)


@MyProject.pre(equilibrated)
@MyProject.post(reacted)
@MyProject.operation
def polymerize(job):
    from scripts.simulate import Simulator
    sinit = Simulator(job)
    sinit.polymerize()
    print("reacted",job.id)

@MyProject.operation
def progress():
    from scripts.progress import progress
    progress()

# Analysis



if __name__ == "__main__":
    MyProject().main()
