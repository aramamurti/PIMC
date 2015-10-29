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

void setup_outfiles(int world_rank, std::ofstream& f1, std::ofstream& f2,std::ofstream& f3,std::ofstream& f4){
    std::stringstream sstm;
    sstm << "overview_" << world_rank <<".txt";
    std::string result = sstm.str();
    f1.open(result.c_str());
    
    std::stringstream sstm2;
    sstm2 << "energydata_" << world_rank <<".csv";
    std::string result2 = sstm2.str();
    f2.open(result2.c_str());
    
    std::stringstream sstm3;
    sstm3 << "permcycledata_" << world_rank <<".csv";
    std::string result3 = sstm3.str();
    f3.open(result3.c_str());
    
    std::stringstream sstm4;
    sstm4 << "windingdata_" << world_rank <<".csv";
    std::string result4 = sstm4.str();
    f4.open(result4.c_str());
}

int main(int argc, const char * argv[]) {
    

    
    MPI_Init(NULL, NULL);
    
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    
    
    std::ofstream f1,f2,f3,f4;
    
    setup_outfiles(world_rank, f1, f2, f3,f4);

    std::cout << world_rank << ": Setting up paths and permutation table..." << std::endl;
    paths* path = new paths(world_rank, f1);
    std::cout << world_rank << ": Started  MC process ..." << std::endl;
    pimc* sim = new pimc();
    
    
    std::vector<double> energy(0);
    std::vector< std::vector<int>> cycles;
    sim->run(path->getParam()->getNumSteps(), path, f2, f3, f4, energy, cycles);
    
    f1 << "Total Energy = " <<path->getUte()->vecavg(energy) << " +/- "<< path->getUte()->vecstd(energy)/sqrt(energy.size())<< std::endl;
    f1 << "Energy/atom= " <<path->getUte()->vecavg(energy)/path->getParam()->getNumParticles()<< "\n" <<std::endl;
    
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
    
    f1 << "Permutation fraction" << std::endl;
    for(std::vector<double>::iterator it = cyclepercent.begin(); it != cyclepercent.end(); it++){
        f1 << it-cyclepercent.begin()+1 <<":\t" << *it <<std::endl;
    }

    f1.close();
        
    delete sim;
    delete path;

    MPI_Finalize();
    
    return 0;
}

