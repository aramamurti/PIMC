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
    sstm2 << "energydata" << world_rank <<".csv";
    std::string result2 = sstm2.str();
    std::ofstream f2;
    f2.open(result2.c_str());
    
    std::stringstream sstm3;
    sstm3 << "permcycledata" << world_rank <<".csv";
    std::string result3 = sstm3.str();
    std::ofstream f3;
    f3.open(result3.c_str());

    paths* path = new paths(world_rank, f);
    std::cout << "Started process " << world_rank << std::endl;
    pimc* sim = new pimc();
    
    
    std::vector<double> energy(0);
    std::vector< std::vector<int>> cycles;
    sim->run(path->getParam()->getNumSteps(), path, f2, f3, energy, cycles);
    
    f << "Total Energy = " <<path->getUte()->vecavg(energy) << " +/- "<< path->getUte()->vecstd(energy)/sqrt(energy.size())<< std::endl;
    f << "Energy/atom= " <<path->getUte()->vecavg(energy)/path->getParam()->getNumParticles()<< "\n" <<std::endl;
    
    int sum = 0;
    for(std::vector<std::vector<int>>::iterator it = cycles.begin(); it != cycles.end(); it++){
        std::vector<int> stepcyc = *it;
        sum += std::accumulate(stepcyc.begin(), stepcyc.end(), 0);
    }
    
    std::vector<int> cyclesum(cycles[0].size(),0);
    for(std::vector<std::vector<int>>::iterator it = cycles.begin(); it != cycles.end(); it++){
        std::vector<int> stepcyc = *it;
        for(std::vector<int>::iterator it2 = stepcyc.begin(); it2 != stepcyc.end(); it2++){
            cyclesum[it2-stepcyc.begin()] += *it2;
        }
    }
    
    std::vector<double> cyclepercent;
    for(std::vector<int>::iterator it = cyclesum.begin(); it != cyclesum.end(); it++){
        cyclepercent.push_back(((double)*it)/sum);
    }
    
    f << "Permutation fraction" << std::endl;
    for(std::vector<double>::iterator it = cyclepercent.begin(); it != cyclepercent.end(); it++){
        f << it-cyclepercent.begin()+1 <<":\t" << *it <<std::endl;
    }

    f.close();
        
    delete sim;
    delete path;

    MPI_Finalize();
    
    return 0;
}

