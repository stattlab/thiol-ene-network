import sys, os, re
import numpy as np
import signac
import gsd, gsd.hoomd 
import networkx as nx 
from collections import defaultdict
import json
from collections import Counter

project = signac.get_project()

for job in project: 
    if job.isfile('polymerize.gsd'):
        # only look at our base cases for now...
        with open(job.fn('signac_statepoint.json')) as f:
            statepoint = json.load(f)
            '''
            if not (statepoint["crosslinker_percent"]==50.0 and statepoint["chain_side_reaction_probability"]==0):
                continue
            '''
            
            if (statepoint["chain_side_reaction_probability"]==0):
                continue
            jobspace = "./workspace/" + str(job) + "/"
            for f in os.listdir(jobspace):
                if not os.path.isfile(jobspace + f):
                    continue 
                if f != "equi.gsd" and f != "signac_job_document.json" and f != "signac_statepoint.json":
                    os.remove(jobspace + f)