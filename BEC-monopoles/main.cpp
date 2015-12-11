//
//  main.cpp
//  PIMC
//
//  Created by Adith Ramamurti on 5/20/15.
//  Copyright (c) 2015 Adith Ramamurti. All rights reserved.
//


/***********************/
 
#include "uni_header.h"
#include "path.h"
#include "mpi.h"
#include "pimc.h"
#include "IO.hpp"

/********************************************************************************************************
 
                                        MAIN METHOD
 
 Method Description: This method runs the simulation. It creates an object to read/write parameters
    and results, sets up a set of paths (world lines) with the given parameters, creates a simulation
    object, and runs the simulation. Finally, the average energy and acceptance ratios are written to a
    file.
 
 ********************************************************************************************************/
int main(int argc, const char * argv[]) {
    
    //Set parameters file
    std::string parameters_file;
    if(argc == 2)
        parameters_file = argv[1];
    else
        parameters_file = "parameters.cfg";
        
    //MPI_Init(NULL, NULL);
    
    int world_rank = 0;
    //MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    
    std::cout << world_rank << ": Setting up..." << std::endl;
    
    //IO object which reads and writes to file
    IO writer(world_rank);

    //Creates path object with input parameters
    boost::shared_ptr<Path> path(new Path(world_rank, writer, writer.read_parameters(parameters_file)));
    
    //Creates a simulation object with the path
    boost::shared_ptr<PIMC> sim(new PIMC(path));
    
    //Initialize energy vector and permutation cycles vector. (Results to be used in analysis)
    dVector energy(0);
    iiVector cycles;
    
    //Run algorithm and get acceptance ratios of moves
    std::cout<< world_rank << ": Starting algorithm..." <<std::endl;
    iVector accept = sim->run(path->get_parameters()->get_end_step(), writer, energy, cycles);
    accept.push_back(path->get_parameters()->get_end_step());
    
    std::cout<< world_rank <<": Finished simulation. Writing results to file..."<< std::endl;
    //Write final results to file
    writer.write_final(path->get_util()->vecavg(energy), path->get_util()->vecstd(energy)/sqrt(energy.size()), path->get_beads()->get_num_particles(),cycles, accept);
    writer.close();
    
    std::cout<< world_rank <<": Done."<< std::endl;
    
    //MPI_Finalize();
    
    return 0;
}

