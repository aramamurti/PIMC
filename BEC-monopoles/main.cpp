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
//#include "mpi.h"
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
    
    int world_rank = 0;
    if(argc == 3){
        std::string parameters_file_dir = argv[1];
        world_rank = std::stoi(argv[2]);
        std::stringstream sstm;
        sstm << parameters_file_dir << "parameters_" << world_rank <<".cfg";
        parameters_file = sstm.str();
    }
    else if(argc == 2)
        parameters_file = argv[1];
    else
        parameters_file = "parameters.cfg";
    
    std::cout << world_rank << ":\tSetting up..." << std::endl;
    
    //IO object which reads and writes to file
    IO writer(world_rank);

    //Creates path object with input parameters
    boost::shared_ptr<Path> path(new Path(world_rank, writer, writer.read_parameters(parameters_file)));
    
    //Creates a simulation object with the path
    boost::shared_ptr<PIMC> sim(new PIMC(path));
    
    //Initialize energy vector and permutation cycles vector. (Results to be used in analysis)
    dVector energy(0);
    iiVector cycles(0);
    dVector particles(0);
    
    //Run algorithm and get acceptance ratios of moves
    std::cout<< world_rank << ":\tStarting algorithm..." <<std::endl;
    std::vector<boost::tuple<std::string, int, int> > accept = sim->run(path->get_parameters()->get_end_step(), writer, energy, cycles, particles);
    
    std::cout<< world_rank <<":\tFinished simulation. Writing results to file..."<< std::endl;
    //Write final results to file
    writer.write_final(path->get_util()->vecavg(energy), path->get_util()->vecstd(energy)/sqrt(energy.size()), path->get_beads()->get_num_particles(),cycles, accept,path->get_util()->vecavg(particles));
    writer.close();
    
    std::cout<< world_rank <<":\tDone."<< std::endl;
    
    return 0;
}

