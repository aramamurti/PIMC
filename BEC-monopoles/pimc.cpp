//
//  PIMC.cpp
//  PIMC
//
//  Created by Adith Ramamurti on 5/20/15.
//  Copyright (c) 2015 Adith Ramamurti. All rights reserved.
//

#include "pimc.h"


PIMC::PIMC(){}

iVector PIMC::run(int end_step, boost::shared_ptr<Path> path, IO &writer, fVector &energytr, iiVector &cycleList){
    
    iVector accept;
    
    std::vector<bool> move_list;
    move_list.push_back(true);
    move_list.push_back(false);
    move_list.push_back(true);
    
    set_up_moves(path, move_list);
    
    for(int step = 0; step < end_step; step++){
        int ptcl = (int) path->get_util()->randnormed(path->get_parameters()->get_num_particles())%path->get_parameters()->get_num_particles();

        moves[0].attempt(ptcl);
        moves[1].attempt(ptcl);
        
        if(step % 100 == 0 && step < path->get_parameters()->get_equilibration()){
            std::cout << path->getPNum() << ": " << step << ", " <<path->energy() << std::endl;
        }
        
        if(step % path->get_parameters()->get_skip() == 0 && step >= path->get_parameters()->get_equilibration()){
            double en = path->energy();
            iVector cycles = path->get_cycles();
            
            energytr.push_back(en);
            cycleList.push_back(cycles);
            
            writer.write_step_state(step, en, path->kinetic_energy(), path->potential_energy(), cycles, path->get_parameters()->get_num_particles(), path->get_winding_number());
        }
    }
    
    for(boost::ptr_vector<Move_Base>::iterator it = moves.begin(); it != moves.end(); it++){
        accept.push_back((*it).get_num_accepts());
    }
    return accept;
    
}

void PIMC::set_up_moves(boost::shared_ptr<Path> path, std::vector<bool> move_list){
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