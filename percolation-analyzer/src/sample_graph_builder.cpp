/*
 * SPDX-FileCopyrightText: 2020 Kevin Höllring for PULS Group <kevin.hoellring@fau.de>
 *
 * SPDX-License-Identifier: MIT
 *
 * Copyright (c) 2020 Kevin Höllring for PULS Group <kevin.hoellring@fau.de>
 *
 * Authors: 2020 Kevin Höllring <kevin.hoellring@fau.de>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the “Software”), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice (including the next paragraph)
 * shall be included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
 * INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
 * PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 */

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <random>
#include <chrono>
#include <fstream>

#include "molecular-graph.hpp"
#include "GSDReadWrite/GSDReader.h"
#include "GSDReadWrite/GSDWriter.h"
#include "GSDReadWrite/VectorMath.h"
#include "GSDReadWrite/SnapshotSystemData.h"
using namespace linalg::aliases;
/**
 * @brief Sample program to randomly generate a molecular graph and convert it into the percolation graph
 *
 * In this program, we will generate a random set of points to populate a cuboid space and generate random bonds between those to build a sample graph to run the percolation analyzer on.
 * This can be replaced by a reading in routine for your respective trajectory format to populate the graph as well as the links/bonds between the atoms by one frame's data.
 * It is mainly meant to provide a sample for how to convert the atomic position and bond data into the percolation graph
 *
 */



 //using namespace linalg::aliases;

 class InputParser{
     public:
         InputParser (int &argc, char **argv){
             for (int i=1; i < argc; ++i)
                 this->tokens.push_back(std::string(argv[i]));
         }

         const std::string& getCmdOption(const std::string &option) const{
             std::vector<std::string>::const_iterator itr;
             itr =  std::find(this->tokens.begin(), this->tokens.end(), option);
             if (itr != this->tokens.end() && ++itr != this->tokens.end()){
                 return *itr;
             }
             static const std::string empty_string("");
             return empty_string;
         }

         bool cmdOptionExists(const std::string &option) const{
             return std::find(this->tokens.begin(), this->tokens.end(), option)
                    != this->tokens.end();
         }
     private:
         std::vector <std::string> tokens;
 };



 /* @param argc
 * @param argv
 * @return int
 */
int main(int argc, char *argv[])
{
    InputParser input(argc, argv);
    if(input.cmdOptionExists("-h")){
        std::cout<<" Percolation algorithm. Usage: sample_graph -i input_gsd -o output_txt" <<std::endl;
        return 0;
    }
    const std::string &filename = input.getCmdOption("-i");
    if (filename.empty()){
        std::cout << "ERROR: need input file, specify '-i input_gsd'"<<std::endl;
        return 1;
    }


    uint64_t frame = 0;
    auto reader = std::make_shared<GSDReader>(filename, frame, false);
    uint64_t Nframes = reader-> getNFrame();

    // 0 Nframes
    int t = 0;
    // Output file
    std::ofstream outputfile;
    outputfile.open(argv[4]);
    // header output
    std::cout<< "# frame percolation_dimension largest_cluster num_molecules" << std::endl;
    outputfile << "# frame percolation_dimension largest_cluster num_molecules" << std::endl;
    for (unsigned int frame_num = 0; frame_num < Nframes; frame_num++)
    {
        auto t_start = std::chrono::high_resolution_clock::now();
        auto reader = std::make_shared<GSDReader>(filename, frame_num, false);
        std::shared_ptr< SnapshotSystemData<float> > snap = reader->getSnapshot();

        //auto ppa = std::make_shared<PrimitivePath>(snap);
        //ppa->removeNonIntersectingTriangles();
        //SnapshotSystemData<float> out_snap = ppa->getSnapshot();
        //writer->analyze(t,out_snap);
        t++;

        //auto t_end = std::chrono::high_resolution_clock::now();
        //double elapsed_time_ms = std::chrono::duration<double, std::milli>(t_end-t_start).count();
        //std::cout << "frame "<< frame_num << "/"<< Nframes <<"  "<< std::endl;

    double3 box = snap->global_box;


    // These are parameters for the generation of the random graph.
    // If you read in your graph, you do not need these
    std::vector< float3 > pos = snap->pos;
    size_t num_points = pos.size();

    double dimx = box[0];
    double dimy = box[1];
    double dimz = box[2];

    double volume = dimx * dimy * dimz;

    //
    // If you read in your molecular graph from a trajectory, this is where to start setting up
    // The respective data structure:
    //
    mol::MolecularGraph mol_graph;

    // Reserve memory for the number of atoms in the frame
    mol_graph.set_atom_count(num_points);

    // Set up a cuboid basis set for the system
    // Can be made triclinic as long as basis[0][1] == basis[0][2] == basis[1][2] == 0
    // as well as basis[0][0] != 0,  basis[1][1] != 0, basis[2][2] != 0
    // The precise setup of this will depend on your chosen basis

    //std::cout << "Register the system pbc" << std::endl;
    std::vector<vec<double>> basis(3);
    basis[0][0] = dimx;
    basis[1][1] = dimy;
    basis[2][2] = dimz;
    //std::cout << "box "<< box[0]<<" " << box[1] << " " << box[2] <<std::endl;
    // Check the return value of set_base to make sure that the graph basis is correctly formatted
    if (!mol_graph.set_basis(basis))
    {
        std::cerr << "Error setting the molecular graph basis. Basis is not triclinic." << std::endl;
        exit(1);
    }


    // Generate (or read in) the position information of the atoms
    std::vector<vec<double>> positions(num_points);

    //
    // If you read in your graph data, this is where you populate the graph with your data instead of randomly generated positions and links.
    //
    for (size_t i = 0; i < num_points; i++)
    {
        positions[i].x = pos[i][0] + box[0]/2.;
        positions[i].y = pos[i][1] + box[1]/2.;
        positions[i].z = pos[i][2] + box[2]/2.;

        mol_graph.set_atom_position(i, positions[i]);

    }
    auto bond_group =  snap->bond_group;
    auto bond_typeid =  snap->bond_type;
    for(int i=0; i<bond_group.size(); i++)
        {
         // don't read dummy bonds
         if(bond_typeid[i]<5)
            {
            mol_graph.add_bond(bond_group[i][0], bond_group[i][1]);
            }
         //std::cout << " bond "<< bond_group[i][0] << " "<< bond_group[i][1] <<std::endl;
     }
    // std::cout << "read "<<  bond_group.size() << " bonds between " << pos.size()<< " particles" << std::endl;
    //
    // At this point, mol_graph holds the entire graph data.
    // We will now convert it to the percolation detection graph
    //

    //std::cout << "Convert molecular graph to percolation graph" << std::endl;
    // Retrieve the Percolation Graph associated with the molecular graph
    percolation::PercolationGraph percolation_graph = mol_graph.get_percolation_graph();

    //std::cout << "Retrieve percolation data" << std::endl;
    // Now calculate the percolation info
    std::vector<percolation::ComponentInfo> percolation_data = percolation_graph.get_component_percolation_info();

    // Let us output the percolation info

    //std::cout << "Results of the percolation analysis:" << std::endl;

//    std::cout << "The system has a total of " << percolation_data.size() << " disjoint molecules" << std::endl;
    //std::cout << "We will now list the data for each of these molecules..." << std::endl
    //          << std::endl;
    bool print = false;
    for (size_t c_index = 0; c_index < percolation_data.size(); c_index++)
    {
        if (percolation_data[c_index].percolation_dim != 0 )
        {

        //std::cout << "Molecule #" << c_index << " stats:" << std::endl;

    //    std::cout << "Percolation dimension: " << percolation_data[c_index].percolation_dim << std::endl;
    //    std::cout << "#Atoms linked: " << percolation_data[c_index].vertices.size() << std::endl;
        std::cout << frame_num << " "<< percolation_data[c_index].percolation_dim  << " "<< percolation_data[c_index].vertices.size() << " " <<percolation_data.size() << std::endl;
        outputfile << frame_num << " "<< percolation_data[c_index].percolation_dim  << " "<< percolation_data[c_index].vertices.size() << " " <<percolation_data.size() << std::endl;
        print = true;
        // If you want full info on which atoms are linked, uncomment the following:
        /*auto &c_verts = percolation_data[c_index].vertices;
        std::cout << "Indices of atoms in this molecule:" << std::endl;
        for (size_t v_index = 0; v_index < c_verts.size(); v_index++)
        {
            std::cout << "\t" << c_verts[v_index].index << std::endl;
        }*/
    }

    }
    if (print==false)
    {
             std::cout << frame_num << " "<< 0 << " " << 0 << " "<< 0 <<std::endl;
             outputfile << frame_num << " "<< 0 << " " << 0 << " "<< 0 <<std::endl;
    }

    } // loop over  frames
}
