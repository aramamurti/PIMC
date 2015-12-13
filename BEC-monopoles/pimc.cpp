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
        
        int num_updates = std::max(1,path->get_beads()->get_num_particles()+path->get_beads()->get_worm_size()/path->get_parameters()->get_num_timeslices());
        for(int i = 0; i < num_updates; i++)
            for(boost::ptr_vector<Move_Base>::iterator it = moves.begin(); it != moves.end(); it++){
                if(it->is_worm_move()){
                    if(!it->is_worm_nec()&&!path->worm_exists())
                        it->attempt();
                    else if(it->is_worm_nec() && path->worm_exists())
                        it->attempt();
                }
                else if(path->get_beads()->get_num_particles() > 0)
                    it->attempt();
            }
        
        if(!path->worm_exists()){
            
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
    
    int off_diag_att = 2000;
    int off_diag_att_ctr = 0;
    
    int conf_att = 200;
    int conf_counter = 0;
    int diag_counter = 0;
    
    int cum_part = 0;
    int start_part = path->get_beads()->get_num_particles();
    
    for(int step = 0; step < end_step; step++){
        std::cout << path->get_processor_num()<<":\t" << step <<"\t"<<estimators[0].estimate()[0] <<std::endl;
        if(double(step)/end_step < 1/4.){
            for(int i = 0; i < std::max(1,path->get_beads()->get_num_particles()+path->get_beads()->get_worm_size()/path->get_parameters()->get_num_timeslices()); i++)
                for(boost::ptr_vector<Move_Base>::iterator it = moves.begin(); it != moves.end(); it++){
                    if(!it->is_worm_move()){
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
            for(int i = 0; i < std::max(1,path->get_beads()->get_num_particles()+path->get_beads()->get_worm_size()/path->get_parameters()->get_num_timeslices()); i++){
                for(boost::ptr_vector<Move_Base>::iterator it = moves.begin(); it != moves.end(); it++){
                    if(it->is_worm_move()){
                        if(!it->is_worm_nec()&&!path->worm_exists())
                            it->attempt();
                        else if(it->is_worm_nec() && path->worm_exists())
                            it->attempt();
                    }
                    else if(path->get_beads()->get_num_particles() > 0)
                        it->attempt();
                }
                
                
                if(double(step)/end_step < 2/3.){
                    cum_part += path->get_beads()->get_num_particles() + path->get_beads()->get_worm_size()/path->get_parameters()->get_num_timeslices();
                    off_diag_att_ctr++;
                    conf_counter++;
                    if(!path->worm_exists())
                        diag_counter++;
                    if(off_diag_att_ctr == off_diag_att){
                        double avg_part = double(cum_part)/off_diag_att_ctr;
                        double mu_shift = 1 - avg_part/start_part;
                        path->get_parameters()->shift_mu(mu_shift);
                        
                        cum_part = 0;
                        off_diag_att_ctr = 0;
                        
                        std::cout << path->get_processor_num() <<":\tChemical pot. equil. -- Avg. particles:\t" << avg_part <<"\t--\tmu:\t" << path->get_parameters()->get_mu() << std::endl;
                        
                    }
                    if(conf_counter == conf_att){
                        double diag_frac = double(diag_counter)/conf_counter;
                        if (path->get_parameters()->get_C0() > 1.0E-5) {
                            if (diag_frac < 0.2)
                                path->get_parameters()->shift_C0(-0.5);
                            else if (diag_frac >= 0.2 && diag_frac < 0.3)
                                path->get_parameters()->shift_C0(-0.4);
                            else if (diag_frac >= 0.3 && diag_frac < 0.4)
                                path->get_parameters()->shift_C0(-0.3);
                            else if (diag_frac >= 0.4 && diag_frac < 0.5)
                                path->get_parameters()->shift_C0(-0.2);
                            else if (diag_frac >= 0.5 && diag_frac < 0.6)
                                path->get_parameters()->shift_C0(-0.1);
                            else if (diag_frac >= 0.6 && diag_frac < 0.75)
                                path->get_parameters()->shift_C0(-0.05);
                        }
                        
                        if (path->get_parameters()->get_C0() < 1.0E4) {
                            if (diag_frac <= 0.9 && diag_frac > 0.85)
                                path->get_parameters()->shift_C0(0.05);
                            else if (diag_frac <= 0.95 && diag_frac > 0.9)
                                path->get_parameters()->shift_C0(0.1);
                            else if (diag_frac > 0.95)
                                path->get_parameters()->shift_C0(0.2);
                        }
                        
                        conf_counter = 0;
                        diag_counter = 0;
                        std::cout << path->get_processor_num()<<":\tWorm equil. -- Diagonal:Total ratio:\t"<< diag_frac <<"\t--\tC0:\t" << path->get_parameters()->get_C0() << std::endl;
                    }
                }
            }
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
                case 3:
                    moves.push_back(new Insert(path));
                    moves.push_back(new Open(path));
                    moves.push_back(new Advance_Head(path));
                    moves.push_back(new Advance_Tail(path));
                    moves.push_back(new Recede_Head(path));
                    moves.push_back(new Recede_Tail(path));
                    moves.push_back(new Swap_Head(path));
                    moves.push_back(new Swap_Tail(path));
                    moves.push_back(new Close(path));
                    moves.push_back(new Remove(path));
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

