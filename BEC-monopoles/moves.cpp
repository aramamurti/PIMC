//
//  Moves.cpp
//  PIMC
//
//  Created by Adith Ramamurti on 5/14/15.
//  Copyright (c) 2015 Adith Ramamurti. All rights reserved.
//

#include "moves.h"


/************************************************************
 
 MOVE BASE CLASS
 
 ************************************************************/

Move_Base::Move_Base(boost::shared_ptr<Path> path){
    num_accepts = 0;
    num_attempts = 0;
    this->path = path;
    pa = boost::shared_ptr<Potential_Action>(new Potential_Action(this->path, this->path->get_parameters()->get_potentials()));
    ndim = path->get_parameters()->get_ndim();
}

void Move_Base::attempt(){
    
    changed_particles.resize(0);
    path->get_beads()->set_old();
    num_particles = path->get_beads()->get_num_particles();
        
    if(num_particles > 0){
        ptcl = (int) path->get_util()->randnormed(num_particles)%num_particles;
    }
    num_attempts++;
    old_action = 0.0;
    new_action = 0.0;
    
}

bool Move_Base::check_move(){
    if(path->get_util()->randnormed(1)<exp(-(new_action - old_action)))
        return true;
    else
        return false;
}

void Move_Base::accept(){
    path->get_beads()->confirm_sep_update();
    num_accepts++;
}

void Move_Base::reject(){
     // std::cout << "Rejected" << std::endl;
    path->get_beads()->revert_old();
}

void Move_Base::reset_acceptance_counters(){
    num_attempts = 0;
    num_accepts = 0;
}

/************************************************************
 
 CENTER OF MASS MOVE CLASS
 
 ************************************************************/

Center_of_Mass::Center_of_Mass(boost::shared_ptr<Path> path) : Move_Base(path){
    
    if(path->get_parameters()->get_box_size() == -1)
        delta = 0.5;
    else{
        delta = 2.*sqrt(path->get_parameters()->get_lambda()*path->get_parameters()->get_tau());
    }
    shift.resize(path->get_parameters()->get_ndim());
    perm_move = false;
    
    move_name = "Center of Mass";
}

void Center_of_Mass::shift_delta(double shift){
    delta += shift*delta;
}

void Center_of_Mass::attempt(){
    Move_Base::attempt();
    if(num_particles == 0){
        reject();
        return;
    }
    
    // std::cout << "Center of Mass Attempt" << std::endl;
    
    
    for(dVector::iterator it = shift.begin(); it != shift.end(); it++ )
        *it = path->get_util()->randgaussian(delta);
    
    changed_particles = path->get_beads()->get_changed_ptcls(ptcl);
    
    for(iVector::iterator it = changed_particles.begin(); it != changed_particles.end(); it++){
        for(int slice = 0; slice < path->get_parameters()->get_num_timeslices(); slice++){
            old_action += pa->get_action_single_particle(*it, slice);
        }
    }
    
    path->get_beads()->shift_all(ptcl, shift, path->get_parameters()->get_box_size());
    
    for(iVector::iterator it = changed_particles.begin(); it != changed_particles.end(); it++){
        for(int slice = 0; slice < path->get_parameters()->get_num_timeslices(); slice++){
            new_action += pa->get_action_single_particle(*it, slice);
        }
    }
    
    if(check_move())
        accept();
    else
        reject();
    
}

void Center_of_Mass::accept(){
    if(path->get_parameters()->get_move_list()[2]){
        for(iVector::iterator it = changed_particles.begin(); it != changed_particles.end(); it++){
            path->push_last_changed(*it);
        }
    }
    Move_Base::accept();
}

/************************************************************
 
 BISECTION MOVE CLASS
 
 ************************************************************/

Bisection::Bisection(boost::shared_ptr<Path> path) : Move_Base(path){
    multistep_dist = path->get_multistep_dist();
    perm_move = false;
    
    move_name = "Bisection";
}

void Bisection::attempt(){
    
    // std::cout << "Bisection Attempt" << std::endl;
    
    Move_Base::attempt();
    
    if(num_particles == 0){
        reject();
        return;
    }
    
    start = path->get_util()->randint(path->get_parameters()->get_num_timeslices());
    
    
    for(int a = 0; a < multistep_dist+1; a++){
        int slice = start+a;
        int row = ptcl;
        if(slice >= path->get_parameters()->get_num_timeslices()){
            slice = slice%path->get_parameters()->get_num_timeslices();
            row = path->get_beads()->get_pair_rows(ptcl, start, multistep_dist)[1];
        }
        old_action += pa->get_action_single_particle(row, slice);
    }
    
    iVector chd_ptcls = path->get_beads()->get_pair_rows(ptcl, start, multistep_dist);
    changed_particles.push_back(chd_ptcls[0]);
    changed_particles.push_back(chd_ptcls[1]);

    level_move(ptcl, start, multistep_dist);
    
    for(int a = 0; a < multistep_dist+1; a++ ){
        int slice = start+a;
        int row = ptcl;
        if(slice >= path->get_parameters()->get_num_timeslices()){
            slice = slice%path->get_parameters()->get_num_timeslices();
            row = path->get_beads()->get_pair_rows(ptcl, start, multistep_dist)[1];
        }
        new_action += pa->get_action_single_particle(row, slice);
    }
    
    if(check_move())
        accept();
    else
        reject();
    
}

void Bisection::level_move(int ptcl, int start, int m){
    if(m != 1 && m%2 == 0){
        int slice = (start + m/2);
        double tau1 = (m/2)*path->get_parameters()->get_tau();
        dVector move;
        move.reserve(ndim);
        dVector bead1;
        dVector bead2;
        path->get_beads()->get_pair_same_path(ptcl, start, m, bead1, bead2);
        
        dVector aved;
        aved.reserve(ndim);
        path->get_util()->avedist(bead1, bead2, aved,path->get_parameters()->get_box_size());
        double width = sqrt(path->get_parameters()->get_lambda()*tau1);
        for(int n = 0; n < ndim; n++){
            move.push_back(aved[n] + path->get_util()->randgaussian(width));
        }
        path->get_beads()->set_bead_data(ptcl, slice, path->get_parameters()->get_box_size(), move);
        level_move(ptcl, start, m/2);
        level_move(ptcl, slice, m/2);
    }
}

void Bisection::accept(){
    std::sort(changed_particles.begin(), changed_particles.end());
    auto last = std::unique(changed_particles.begin(), changed_particles.end());
    changed_particles.erase(last, changed_particles.end());
    
    if(path->get_parameters()->get_move_list()[2]){
        for(iVector::iterator it = changed_particles.begin(); it != changed_particles.end(); it++){
            path->push_last_changed(*it);
        }
        path->set_last_start_end(start, start+multistep_dist);
        path->get_beads()->set_prev_perm();
    }
    Move_Base::accept();
}

void Bisection::reject(){
    Move_Base::reject();
}

/**********************************************************************************
 
 PERMUTATION BISECTION MOVE CLASS
 
 ***********************************************************************************/

Perm_Bisection::Perm_Bisection(boost::shared_ptr<Path> path) : Bisection(path){
    ptable = boost::shared_ptr<Permutation_Table>(new Permutation_Table(path));
    perm_move = true;
    
    move_name = "Permutation Bisection";
    
}

void Perm_Bisection::attempt(){
    Move_Base::attempt();
    
    if(num_particles == 0){
        reject();
        return;
    }
    
    
    start = path->get_util()->randint(path->get_parameters()->get_num_timeslices());
    
    for(int slice = 0; slice < path->get_parameters()->get_num_timeslices(); slice++)
        ptable->recalc_perms(path->get_last_changed(), slice);
    
    
    for(int a = 0; a < multistep_dist+1; a++){
        int slice = (start + a)%path->get_parameters()->get_num_timeslices();
        old_action += pa->get_action(slice,0);
    }
    
    iVector identity(num_particles);
    iota(identity.begin(),identity.end(),0);
    iVector chosenPerm = identity;
    iVector origpart(0);
    permed_parts = std::vector<int>(0);
    
    if(path->get_parameters()->is_boson()){
        ptable->pick_permutation(start, chosenPerm);
        
        for(iVector::iterator it = identity.begin(); it != identity.end(); it++)
            if(*it != chosenPerm[*it]){
                origpart.push_back(*it);
                permed_parts.push_back(chosenPerm[*it]);
            }
    }
    
    if(permed_parts.size() == 0){
        iVector chd_ptcls = path->get_beads()->get_pair_rows(ptcl, start, multistep_dist);
        changed_particles.push_back(chd_ptcls[0]);
        changed_particles.push_back(chd_ptcls[1]);

        level_move(ptcl, start, multistep_dist);
        permed_parts.push_back(ptcl);
    }
    else{
        path->get_beads()->set_permutation(origpart, permed_parts, start, multistep_dist);
        path->get_beads()->permute();
        for(iVector::iterator ptcl = origpart.begin(); ptcl !=origpart.end(); ptcl++){
            iVector chd_ptcls = path->get_beads()->get_pair_rows(*ptcl, start, multistep_dist);
            changed_particles.push_back(chd_ptcls[0]);
            changed_particles.push_back(chd_ptcls[1]);

            level_move(*ptcl, start, multistep_dist);
        }
    }
    
    for(int a = 0; a < multistep_dist+1; a++ ){
        int slice = (start + a)%path->get_parameters()->get_num_timeslices();
        new_action += pa->get_action(slice,0);
    }
    
    if(check_move())
        accept();
    else
        reject();
}

void Perm_Bisection::accept(){
    
    auto last = std::unique(changed_particles.begin(), changed_particles.end());
    changed_particles.erase(last, changed_particles.end());
    path->set_last_changed(changed_particles);
    path->set_last_start_end(start, start + multistep_dist);
    path->get_beads()->set_prev_perm();
    path->get_beads()->reset_indices();
    Move_Base::accept();
    
}
void Perm_Bisection::reject(){
    if(permed_parts.size() != 0){
        path->get_beads()->permute(true);
        path->get_beads()->reset_permute();
    }
    path->set_last_changed(iVector(0));
    Bisection::reject();
}
