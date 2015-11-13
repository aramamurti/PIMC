//
//  PIMC.cpp
//  PIMC
//
//  Created by Adith Ramamurti on 5/20/15.
//  Copyright (c) 2015 Adith Ramamurti. All rights reserved.
//

#include "pimc.h"


PIMC::PIMC(boost::shared_ptr<Path> path){this->path = path;}

iVector PIMC::run(int end_step, IO &writer, fVector &energytr, iiVector &cycleList){
    
    iVector accept;
    
    std::vector<bool> move_list;
    move_list.push_back(false);
    move_list.push_back(false);
    move_list.push_back(true);
    
    set_up_moves(move_list);
    
    en = boost::shared_ptr<Energy_Estimator>(new Energy_Estimator(path, path->get_parameters()->get_potentials()));
    
    for(int step = 0; step < end_step; step++){
        int ptcl = (int) path->get_util()->randnormed(path->get_parameters()->get_num_particles())%path->get_parameters()->get_num_particles();

        
        if(step % 100 == 0 && step < path->get_parameters()->get_equilibration()){
            std::cout << path->getPNum() << ": " << step << ", " <<en->total_energy() << std::endl;
        }
        moves[0].attempt(ptcl);
        //moves[1].attempt(ptcl);

        
        if(step % path->get_parameters()->get_skip() == 0 && step >= path->get_parameters()->get_equilibration()){
            iVector cycles = path->get_cycles();
            
            energytr.push_back(en->total_energy());
            cycleList.push_back(cycles);
            
            writer.write_step_state(step, en->total_energy(), en->kinetic_energy(), en->potential_energy(), cycles, path->get_beads()->get_num_particles(), path->get_winding_number());
        }
    }
    
    for(boost::ptr_vector<Move_Base>::iterator it = moves.begin(); it != moves.end(); it++){
        accept.push_back((*it).get_num_accepts());
    }
    return accept;
    
}

void PIMC::set_up_moves(std::vector<bool> move_list){
    int i = 0;
    
    for(std::vector<bool>::iterator it = move_list.begin(); it != move_list.end(); it++){
        if(*it)
            switch(i){
                case 0:
                    moves.push_back(new Center_of_Mass(path));
                    break;
                case 1:
                    moves.push_back(new Bisection(path));
                    break;
                case 2:
                    moves.push_back(new Perm_Bisection(path));
                    break;
            }
        i++;
    }
}