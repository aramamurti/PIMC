//
//  PIMC.cpp
//  PIMC
//
//  Created by Adith Ramamurti on 5/20/15.
//  Copyright (c) 2015 Adith Ramamurti. All rights reserved.
//

#include "pimc.h"


PIMC::PIMC(boost::shared_ptr<Path> path){this->path = path;}

std::vector<boost::tuple<std::string, int, int> > PIMC::run(int end_step, IO &writer, dVector &energytr, iiVector &cycleList, dVector &particles){
    
    std::vector<boost::tuple<std::string, int, int> > accept;
    
    std::vector<bool> estimator_list(3,true);
    
    set_up_moves(path->get_parameters()->get_move_list());
    set_up_estimators(estimator_list);
    
    std::cout<< path->get_processor_num() << ":\tStarting algorithm..." <<std::endl;

    std::cout << path->get_processor_num() <<":\tEquilibrating..." <<std::endl;
    
    equilibrate();
    
    for(boost::ptr_vector<Move_Base>::iterator it = moves.begin(); it != moves.end(); it++){
        it->reset_acceptance_counters();
        if(it->get_move_name() == "Center of Mass"){
            writer.write_equil_parameters(path->get_parameters(), it->get_delta());
        }
    }
    
    std::cout << path->get_processor_num() <<":\tStarting simulation..." << std::endl;
    for(int step = 0; step < end_step; step++){
        
        std::cout << path->get_processor_num()<<":\tSimulation Step:\t" << step <<"/"<<end_step <<std::endl;

        int num_updates = std::max(1,path->get_beads()->get_num_particles());
        for(int i = 0; i < num_updates; i++)
            for(boost::ptr_vector<Move_Base>::iterator it = moves.begin(); it != moves.end(); it++){
                if(path->get_beads()->get_num_particles() > 0)
                    it->attempt();
            }
        
        if(step% 20 == 0){
            
            dVector cycles_double = estimators[1].estimate();
            iVector cycles(cycles_double.begin(),cycles_double.end());
            
            dVector winding_double = estimators[2].estimate();
            iVector winding(winding_double.begin(), winding_double.end());
            dVector energy = estimators[0].estimate();
            
            energytr.push_back(energy[0]);
            cycleList.push_back(cycles);
            particles.push_back(path->get_beads()->get_num_particles());
            
            writer.write_step_state(step, energy, cycles, path->get_beads()->get_num_particles(), winding);
        }
    }
    
    for(boost::ptr_vector<Move_Base>::iterator it = moves.begin(); it != moves.end(); it++){
        accept.push_back(boost::tuple<std::string, int, int>(it->get_move_name(), it->get_num_accepts(), it->get_num_attempts()));
    }
    return accept;
    
}

void PIMC::equilibrate(){
    int end_step = path->get_parameters()->get_equilibration();
    
    int com_att = 200;
    int com_acc = 0;
    
    for(int step = 0; step < end_step; step++){
        
        std::cout << path->get_processor_num()<<":\tEquilibration Step:\t" << step <<"/"<<end_step <<std::endl;
        
        if(double(step)/end_step < 1/3.){
            for(int i = 0; i < std::max(1,path->get_beads()->get_num_particles()); i++)
                for(boost::ptr_vector<Move_Base>::iterator it = moves.begin(); it != moves.end(); it++){
                    if(!it->is_perm_move()){
                        it->attempt();
                        if((it->get_move_name() == "Center of Mass") && (it->get_num_attempts()%com_att == 0)){
                            if(!(it->get_delta() > path->get_parameters()->get_box_size()/2.)){
                                com_acc = it->get_num_accepts() - com_acc;
                                double com_acc_rat = double(com_acc)/com_att;
                                if (com_acc_rat < 0.2)
                                    it->shift_delta(-0.6);
                                else if (com_acc_rat < 0.3)
                                    it->shift_delta(-0.4);
                                else if (com_acc_rat < 0.4)
                                    it->shift_delta(-0.2);
                                else if (com_acc_rat > 0.6)
                                    it->shift_delta(0.2);
                                else if (com_acc_rat > 0.7)
                                    it->shift_delta(0.4);
                                else if (com_acc_rat > 0.8)
                                    it->shift_delta(0.6);
                                com_acc = it->get_num_accepts();
                                std::cout << path->get_processor_num() <<":\tCenter of Mass equil. -- Acceptance:\t"<<com_acc_rat <<"\t--\tDelta:\t"<<it->get_delta() << std::endl;
                            }
                        }
                    }
                }
        }
        else{
            for(int i = 0; i < std::max(1,path->get_beads()->get_num_particles()); i++)
                for(boost::ptr_vector<Move_Base>::iterator it = moves.begin(); it != moves.end(); it++)
                        it->attempt();
            
        }
    }
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

