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
    
    MPI_Init(NULL, NULL);
    
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    
    paths* path = new paths(world_rank+1);
    pimc sim;
    vector<bool> pmv = path->getParam()->getMoves();
    vector<double> energy = sim.run(path->getParam()->getNumSteps(), path, pmv);
    cout<< "Energy = " <<path->getUte()->vecavg(energy) << " +/- "<< path->getUte()->vecstd(energy)/sqrt(energy.size())<<"\n";
    delete path;

    MPI_Finalize();
    
    return 0;
}

