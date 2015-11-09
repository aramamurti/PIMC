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
    
    MPI_Init(NULL, NULL);
    
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    
    
    std::ofstream f1,f2,f3,f4;
    
    IO writer(world_rank);

    Path* path = new Path(world_rank, writer);
    Pimc* sim = new Pimc();
    
    vectorf energy(0);
    vectorii cycles;
    
    vectori accept = sim->run(path->get_parameters()->get_end_step(), path, writer, energy, cycles);
    accept.push_back(path->get_parameters()->get_end_step());
    writer.write_final(path->getUte()->vecavg(energy), path->getUte()->vecstd(energy)/sqrt(energy.size()), path->get_parameters()->get_num_particles(),cycles, accept);
    
    delete sim;
    delete path;

    writer.close();
    
    MPI_Finalize();
    
    return 0;
}

