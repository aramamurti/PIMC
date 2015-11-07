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

void write_outfiles(std::ofstream &f1, Path* path, vectorf energy, vectorii cycles){
    
    f1 << "Total Energy = " <<path->getUte()->vecavg(energy) << " +/- "<< path->getUte()->vecstd(energy)/sqrt(energy.size())<< std::endl;
    f1 << "Energy/atom= " <<path->getUte()->vecavg(energy)/path->get_parameters()->get_num_particles()<< "\n" <<std::endl;
    
    int sum = 0;
    for(vectorii::iterator it = cycles.begin(); it != cycles.end(); it++){
        vectori stepcyc = *it;
        sum += std::accumulate(stepcyc.begin(), stepcyc.end(), 0);
    }
    
    vectori cyclesum(cycles[0].size(),0);
    for(vectorii::iterator it = cycles.begin(); it != cycles.end(); it++){
        vectori stepcyc = *it;
        for(vectori::iterator it2 = stepcyc.begin(); it2 != stepcyc.end(); it2++){
            cyclesum[it2-stepcyc.begin()] += *it2;
        }
    }
    
    vectorf cyclepercent;
    for(vectori::iterator it = cyclesum.begin(); it != cyclesum.end(); it++){
        cyclepercent.push_back(((double)*it)/sum);
    }
    
    f1 << "Permutation fraction" << std::endl;
    for(vectorf::iterator it = cyclepercent.begin(); it != cyclepercent.end(); it++){
        f1 << it-cyclepercent.begin()+1 <<":\t" << *it <<std::endl;
    }
    
    f1.close();
}

int main(int argc, const char * argv[]) {
    
    MPI_Init(NULL, NULL);
    
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    
    
    std::ofstream f1,f2,f3,f4;
    
    setup_outfiles(world_rank, f1, f2, f3,f4);

    Path* path = new Path(world_rank, f1);
    Pimc* sim = new Pimc();
    
    vectorf energy(0);
    vectorii cycles;
    
    sim->run(path->get_parameters()->get_end_step(), path, f2,f3,f4, energy, cycles);
    write_outfiles(f1,path,energy,cycles);
    
    delete sim;
    delete path;

    MPI_Finalize();
    
    return 0;
}

