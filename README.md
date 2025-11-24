# thiol-ene-network

Software requirements  
-----------------
modules:  anaconda/3 
conda: python 3.10, numpy, scipy, freud, networkx, omnia::eigen3, pybind11, hoomd 4.0.1, gsd 2.X, signac 2.X, signac-flow

also needed: matplotib etc. the usual
for pretty graph visualization: netgraph grst::rectangle-packer


Run instructions
---------------- 

to initialize or to add more statepoints:
python init.py

to check on progress:
python project.py status

to run locally 
python project.py run -o equilibrate (-n 1)
python project.py run -o polymerize -n 1 -j 7f0fd52c32da20cd80a98cc394647ab8 