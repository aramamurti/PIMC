//
//  pimc.cpp
//  PIMC
//
//  Created by Adith Ramamurti on 5/20/15.
//  Copyright (c) 2015 Adith Ramamurti. All rights reserved.
//

#include "moves.h"
#include "paths.h"
#include "pimc.h"
#include "uni_header.h"
#include <unistd.h>



Pimc::Pimc(){
    numacceptc = 0;
    numaccepts = 0;
    numacceptb = 0;
}

void Pimc::run(int end_step, Path* path, std::ofstream &f1,std::ofstream &f2,std::ofstream &f3, vectorf &energytr, vectorii &cycleList){
    
    
    moves mvs;
        
    f1 << "Step No." << ", " << "Energy/atom" << ", " << "KE/atom" << ", " << "PE/atom" << std::endl;

    
    for(int step = 0; step < end_step; step++){
        int ptcl = (int) path->getUte()->randnormed(path->get_parameters()->get_num_particles())%path->get_parameters()->get_num_particles();
        if(mvs.comMove(path, ptcl))
            numacceptc += 1;
        if(mvs.bisectionMoveHelper(path, ptcl))
            numacceptb += 1;
        
        if(step == 500 && step < path->get_parameters()->get_equilibration()){
            std::cout << path->getPNum() << ": " << step << ", " <<path->energy() << std::endl;
            path->get_beads()->print_list_file(step);

        }
        if(step % path->get_parameters()->get_skip() == 0 && step >= path->get_parameters()->get_equilibration()){
            double en = path->energy();
            f1 << step << ", " << en/path->get_parameters()->get_num_particles() << ", " << path->kineticEnergy()/path->get_parameters()->get_num_particles() << ", " << path->potentialEnergy()/path->get_parameters()->get_num_particles() << std::endl;
            std::vector<int> cycles = path->get_cycles();
            for(std::vector<int>::iterator it = cycles.begin(); it != cycles.end(); it++){
                f2 << *it;
                if(cycles.size() - (it-cycles.begin()) != 1)
                    f2 << ", ";
            }
            f2 << std::endl;
            
            std::vector<int> wnum = path->get_winding_number();
            for(std::vector<int>::iterator it = wnum.begin(); it != wnum.end(); it++){
                f3 << *it;
                if(wnum.size() - (it-wnum.begin()) != 1)
                    f3 << ", ";
            }
            f3 << std::endl;
            
            energytr.push_back(en);
            cycleList.push_back(cycles);
        }
    }
    
    f1 << "\nCenter of mass acceptance: "<< 1.0*numacceptc/(end_step) << "\nBisection acceptance: "<< 1.0*numacceptb/(end_step) << "\n" << std::endl;
    
}