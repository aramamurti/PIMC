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


int main(int argc, const char * argv[]) {
    
    std::string parameters_file;
    if(argc == 2)
        parameters_file = argv[1];
    else
        parameters_file = "parameters.cfg";
        
    MPI_Init(NULL, NULL);
    
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    
    std::cout << world_rank << ": Setting up..." << std::endl;
    
    IO writer(world_rank);

    boost::shared_ptr<Path> path(new Path(world_rank, writer, writer.read_parameters(parameters_file)));
    
    boost::shared_ptr<PIMC> sim(new PIMC(path));
    
    dVector energy(0);
    iiVector cycles;
    
    std::cout<< world_rank << ": Starting algorithm..." <<std::endl;
    iVector accept = sim->run(path->get_parameters()->get_end_step(), writer, energy, cycles);
    accept.push_back(path->get_parameters()->get_end_step());
    writer.write_final(path->get_util()->vecavg(energy), path->get_util()->vecstd(energy)/sqrt(energy.size()), path->get_beads()->get_num_particles(),cycles, accept);

    writer.close();
    
    MPI_Finalize();
    
    return 0;
}

