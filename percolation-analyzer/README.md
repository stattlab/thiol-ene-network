# Percolation Detection Analysis

## Description

In a finite MD system with periodic boundary conditions (pbc) of any (triclinic) shape, it is relevant to detect the percolation point in order to formally define the point of gelation.
Here the percolation point is defined as the time, at which a molecule first percolates (= connects to itself) in each dimension.

This repository contains the necessary code to analyze whether or not a molecule is percolating, in how many dimensions it does so (1d -> string like, 2d -> sheet like, 3d -> grid-like) and returns the information for each individual connected molecule in the graph provided.

This code is build exisiting code from <https://github.com/puls-group/percolation-analyzer.git>
If you use the code provided in this repository for research purposes, please cite the paper "An Exact Algorithm to Detect the Percolation Transition in Molecular Dynamics Simulations of Cross-Linking Polymer Networks" (doi: [10.1021/acs.jctc.1c00423](https://doi.org/10.1021/acs.jctc.1c00423 )) as a reference.

## License/Permission of use

The code is provided as-is without any warranty for its correctness under the MIT license (see `LICENSE` file).

## How to use

`src/sample_graph_builder.cpp` is the main, where the gsd is read and iterated over. The output
goes straight to the terminal and tells you, for each frame:

1. percolation dimension (0-not percolated, 1-1D percolation, 2-2D percolation, 3-3D percolation)
2. size of the percolating cluster
3. number of disconnected clusters in the system

If you want to know all the monomers in the percolating clusters, comment out the code at the end of
`src/sample_graph_builder.cpp`  and recompile.

### Building the library/Build system

You need cmake and a C++/C compiler.

Execute the following commands in the project folder (`percolation-analyzer`):

```
mkdir build
cd build
cmake ..
make
```

This will create an exectuable `bin/sample_graph` that you can use.

## Authors

The code was written by Kevin Höllring in 2020, phd student at PULS group of
Friedrich-Alexander-University Erlangen-Nuremberg, modified by Antonia Statt at UIUC to
work with hoomd-blue GSD file trajectories.
