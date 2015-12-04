//
//  PIMC.cpp
//  PIMC
//
//  Created by Adith Ramamurti on 5/20/15.
//  Copyright (c) 2015 Adith Ramamurti. All rights reserved.
//

#include "pimc.h"


PIMC::PIMC(boost::shared_ptr<Path> path){this->path = path;}

iVector PIMC::run(int end_step, IO &writer, dVector &energytr, iiVector &cycleList){
    
    iVector accept;
    
    std::vector<bool> estimator_list(3,true);
    
    set_up_moves(path->get_parameters()->get_move_list());
    set_up_estimators(estimator_list);
    
    for(int step = 0; step < end_step; step++){
        
        if(step % 100 == 0 && step < path->get_parameters()->get_equilibration()){
            std::cout << path->get_processor_num() << ": " << step << ", " <<(estimators[0].estimate())[0] << std::endl;
        }
        
        for(boost::ptr_vector<Move_Base>::iterator it = moves.begin(); it != moves.end(); it++){
            it->attempt();
        }
        
        if(step % path->get_parameters()->get_skip() == 0 && step >= path->get_parameters()->get_equilibration()){
            dVector cycles_double = estimators[1].estimate();
            iVector cycles(cycles_double.begin(),cycles_double.end());
            
            dVector winding_double = estimators[2].estimate();
            iVector winding(winding_double.begin(), winding_double.end());
            dVector energy = estimators[0].estimate();
            
            energytr.push_back(energy[0]);
            cycleList.push_back(cycles);
            
            writer.write_step_state(step, energy, cycles, path->get_beads()->get_num_particles(), winding);
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
                case 3:
                    moves.push_back(new Insert(path));
                    moves.push_back(new Remove(path));
                    moves.push_back(new Open(path));
                    moves.push_back(new Close(path));
                    moves.push_back(new Advance_Head(path));
                    moves.push_back(new Advance_Tail(path));
                    moves.push_back(new Recede_Head(path));
                    moves.push_back(new Recede_Tail(path));
                    moves.push_back(new Swap_Head(path));
                    moves.push_back(new Swap_Tail(path));
                    break;
            }
        i++;
    }
}

void PIMC::set_up_estimators(std::vector<bool> estimator_list){
    int i = 0;
    for(std::vector<bool>::iterator it = estimator_list.begin(); it != estimator_list.end(); it++){
        if(*it)
            switch(i){
                case 0:
                    estimators.push_back(new Energy_Estimator(path,path->get_parameters()->get_potentials()));
                    break;
                case 1:
                    estimators.push_back(new Permutation_Estimator(path));
                    break;
                case 2:
                    estimators.push_back(new Winding_Estimator(path));
                    break;
            }
        i++;
    }
}

