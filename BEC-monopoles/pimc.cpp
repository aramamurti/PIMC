//
//  pimc.cpp
//  PIMC
//
//  Created by Adith Ramamurti on 5/20/15.
//  Copyright (c) 2015 Adith Ramamurti. All rights reserved.
//

#include "pimc.h"


Pimc::Pimc(){
    numacceptc = 0;
    numaccepts = 0;
    numacceptb = 0;
}

std::vector<int> Pimc::run(int end_step, Path* path, IO &writer, vectorf &energytr, vectorii &cycleList){
    
    
    moves mvs;
        
    std::vector<int> accept;
    
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
            std::vector<int> cycles = path->get_cycles();
            
            energytr.push_back(en);
            cycleList.push_back(cycles);
            
            writer.write_step_state(step, en, path->kineticEnergy(), path->potentialEnergy(), cycles, path->get_parameters()->get_num_particles(), path->get_winding_number());
        }
    }
    
    accept.push_back(numacceptc);
    accept.push_back(numacceptb);
    return accept;
    
}