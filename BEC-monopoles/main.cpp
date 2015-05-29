//
//  main.cpp
//  PIMC
//
//  Created by Adith Ramamurti on 5/20/15.
//  Copyright (c) 2015 Adith Ramamurti. All rights reserved.
//


/***********************/
 
#include <iostream>
#include <vector>
#include "moves.h"
#include "paths.h"
#include "pimc.h"
#include "mpi.h"

using namespace std;


int main(int argc, const char * argv[]) {
    
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);
    
    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    
    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);
    
    // Print off a hello world message
    printf("Running from processor %d"
           " out of %d processors\n", world_rank+1, world_size);
    
    paths* path = new paths(world_rank+1);
    pimc sim;
    vector<bool> pmv = path->getParam()->getMoves();
    vector<double> energy = sim.run(path->getParam()->getNumSteps(), path, pmv);
    cout<< "Energy = " <<path->getUte()->vecavg(energy) << " +/- "<< path->getUte()->vecstd(energy)/sqrt(energy.size())<<"\n";
    delete path;

    
    // Finalize the MPI environment.
    MPI_Finalize();
    
    return 0;
}

