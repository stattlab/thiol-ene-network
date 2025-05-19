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
    cutoff = 0
    if float(conversion_data[-1,1]) - float(conversion_data[-5,1]) <= cutoff:
        return True
    else:
        return False

@MyProject.label
def reacting(job):
    return job.isfile("polymerize.gsd")

@MyProject.label
def deformed(job):
    return job.isfile("deform.gsd") and os.path.getsize(job.fn('deform.log')) > 100
    
@MyProject.label
def contracted_bonds(job):
    return job.isfile("contract_bonds.gsd")

#-----------------------
# Analysis labels
#-----------------------

# def percolation_analyzed(job):
#     return job.isfile('percolation_data.txt')
@MyProject.label
def defect_analyzed(job):
    return job.isfile('loop_counts.txt')

@MyProject.label
def xlink_rdf_analyzed(job):
    #If the system cannot form any crosslinks, don't even try to make an RDF of them
    if job.sp['chain_side_reaction_probability'] != 0 :
        one = job.isfile('ene_ene_rdf.txt')
    else:
        one = True
    if job.sp['crosslinker_percent'] != 0:
        two = job.isfile('thiol_thiol_rdf.txt')
    else:
        two = True
    if job.sp['crosslinker_percent'] != 0 and job.sp['chain_side_reaction_probability'] != 0:
        three = job.isfile('thiol_ene_rdf.txt')
    else:
        three = True
    return one and two and three

@MyProject.label
def stress_strain_analyzed(job):
    return job.isfile('stress_strain.txt')

@MyProject.label
def modulus_analyzed(job):
    return job.isfile('youngs_modulus.txt')

@MyProject.label
def gelation_conversion_2(job):
    try:
        return job.doc['gelation_conversion_2_2'][0] >= 0.0
    except KeyError:
        return False
    
@MyProject.label
def contract_bonds_analyzed(job):
    return job.isfile('contract_bonds_analysis/effective_strand_hist.txt') and \
    job.isfile('contract_bonds_analysis/ineffective_strand_hist.txt') and \
    job.isfile('contract_bonds_analysis/overall.txt') and \
    job.isfile('contract_bonds_analysis/scanlan_case_analysis.txt') and \
    job.isfile('contract_bonds_analysis/crosslink_properties.txt')

@MyProject.label
def vv_analyzed(job):
    return job.isfile('contract_bonds_analysis/voronoi_volumes_all.txt') and \
            job.isfile('contract_bonds_analysis/voronoi_volumes_thiol.txt') and \
            job.isfile('contract_bonds_analysis/voronoi_volumes_ene.txt')



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

@MyProject.pre(reacted)
@MyProject.post(deformed)
@MyProject.operation
def deform(job):
    import scripts.deform
    scripts.deform.main(job.fn('polymerize.gsd'),-1)
    print('deformed: ',job.id)

@MyProject.pre(reacted)
@MyProject.post(contracted_bonds)
@MyProject.operation
def contract_bonds(job):
    from scripts.simulate import Simulator
    sinit = Simulator(job)
    sinit.contract_bonds()
    print("contracted bonds: ",job.id)

#-----------------------
# Analysis operations
#-----------------------
# @MyProject.pre(reacting)
# @MyProject.post(trajectory_conversion_analyzed)
# @MyProject.operation
# def analyze_trajectory_conversion(job):
#     os.system('python3 ./scripts/trajectory_conversion.py')

@MyProject.pre(reacting)
@MyProject.post(reacted)
@MyProject.operation
def calculate_conversion(job):
    import scripts.conversion
    print("now conducting trajectory conversion analysis")
    scripts.conversion.calculate_conversion(job.id)
    print('analyzed the trajectory conversion of ', job.id)

@MyProject.pre(reacted)
@MyProject.post(defect_analyzed)
@MyProject.operation
def defect_analysis(job):
    import scripts.network_properties
    scripts.network_properties.defect_analysis(job.id)
    print('analyzed defectivity: ',job.id)

@MyProject.pre(reacted)
@MyProject.post(xlink_rdf_analyzed)
@MyProject.operation
def xlink_rdf_analysis(job):
    import scripts.network_properties
    scripts.network_properties.xlink_rdf_analysis(job.id)
    print('analyzed crosslink rdfs: ',job.id)


@MyProject.pre(reacted)
@MyProject.operation
def strand_analysis(job):
    import scripts.network_properties
    scripts.network_properties.strand_lengths_analysis(job.id)
    print('analyzed strands: ',job.id)


@MyProject.pre(deformed)
@MyProject.post(stress_strain_analyzed)
@MyProject.operation
def deformation_analysis(job):
    import scripts.plot_deforms
    scripts.plot_deforms.analyze_stress_strain(job.id)
    print('analyzed stress-strain: ',job.id)

@MyProject.pre(stress_strain_analyzed)
@MyProject.post(modulus_analyzed)
@MyProject.operation
def modulus_analysis(job):
    import scripts.plot_deforms
    scripts.plot_deforms.analyze_modulus(job.id)
    print('analyzed modulus: ',job.id)

@MyProject.pre(reacted)
@MyProject.post(gelation_conversion_2)
@MyProject.operation
def gelation_2_analysis(job):
    import scripts.gelation2
    scripts.gelation2.gelation_analysis_via_reducedMw(job.id)
    print('analyzed gelation via reduced Mw method: ',job.id)

@MyProject.pre(contracted_bonds)
@MyProject.post(contract_bonds_analyzed)
@MyProject.operation
def contract_bonds_analysis(job):
    import scripts.network_properties
    scripts.network_properties.contract_bonds_analysis(job.id)
    scripts.network_properties.scanlan_case_analysis_on_contract_bonds(job.id)
    print('analyzed contract bonds: ',job.id)

@MyProject.pre(reacted and contract_bonds_analyzed)
@MyProject.post(vv_analyzed)
@MyProject.operation
def vv_analysis(job):
    import scripts.network_properties
    scripts.network_properties.crosslinker_heterogeneity_by_VV(job.id)
    print('analyzed voronoi volumes of job: ',job.id)

if __name__ == "__main__":
    MyProject().main()
