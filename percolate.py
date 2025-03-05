import os 
import signac 

# our signac project
project = signac.get_project()
for job in project: 
    if job.sp["replica_index"]!=1:
        continue
    polymerize_file = "./workspace/" + str(job.id) + "/polymerize.gsd"
    os.system("./percolation-analyzer/bin/sample_graph -i " + polymerize_file + " -o " + "./workspace/" + str(job.id) + "/percolation_data.txt")