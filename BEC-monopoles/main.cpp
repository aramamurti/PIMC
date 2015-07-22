//
//  main.cpp
//  PIMC
//
//  Created by Adith Ramamurti on 5/20/15.
//  Copyright (c) 2015 Adith Ramamurti. All rights reserved.
//


/***********************/
 
#include "uni_header.h"
#include "paths.h"
#include "mpi.h"
#include "pimc.h"


int main(int argc, const char * argv[]) {
    

    
    MPI_Init(NULL, NULL);
    
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    
    std::stringstream sstm;
    sstm << "overview" << world_rank <<".txt";
    std::string result = sstm.str();
    std::ofstream f;
    f.open(result.c_str());
    
    std::stringstream sstm2;
    sstm2 << "data" << world_rank <<".csv";
    std::string result2 = sstm.str();
    std::ofstream f2;
    f2.open(result.c_str());

    paths* path = new paths(world_rank, f);
    std::cout << "Started process " << world_rank << std::endl;
    pimc* sim = new pimc();
    
    
    
    std::vector<double> energy = sim->run(path->getParam()->getNumSteps(), path, f2);
    f << "Total Energy = " <<path->getUte()->vecavg(energy) << " +/- "<< path->getUte()->vecstd(energy)/sqrt(energy.size())<< std::endl;
    f << "Energy/atom= " <<path->getUte()->vecavg(energy)/path->getParam()->getNumParticles()<< std::endl;

    f.close();
    
    std::cout << path->numswap;
    
    delete sim;
    delete path;

    MPI_Finalize();
    
    return 0;
}

