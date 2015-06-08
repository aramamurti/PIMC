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
#include "mpi.h"
#include "pimc.h"


int main(int argc, const char * argv[]) {
    
    MPI_Init(NULL, NULL);
    
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    
    paths* path = new paths(world_rank+1);
    std::cout << "Started process " << world_rank+1 << std::endl;
    pimc sim;
    std::vector<bool> pmv = path->getParam()->getMoves();
    std::vector<double> energy = sim.run(path->getParam()->getNumSteps(), path, pmv);
    std::cout<< "Energy = " <<path->getUte()->vecavg(energy) << " +/- "<< path->getUte()->vecstd(energy)/sqrt(energy.size())<<"\n";
    delete path;

    MPI_Finalize();
    
    return 0;
}

