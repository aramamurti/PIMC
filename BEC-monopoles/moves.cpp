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
}

void Move_Base::attempt(){
    ptcl = (int) path->get_util()->randnormed(path->get_parameters()->get_num_particles())%path->get_parameters()->get_num_particles();
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
    num_accepts++;
    path->put_in_box();
}

void Move_Base::reject(){
    path->get_beads()->revert_old_data();
}

/************************************************************
 
                CENTER OF MASS MOVE CLASS
 
 ************************************************************/

Center_of_Mass::Center_of_Mass(boost::shared_ptr<Path> path) : Move_Base(path){
    
    if(path->get_parameters()->get_box_size() == -1)
        delta = 0.5;
    else{
        delta = sqrt(path->get_parameters()->get_lambda()*path->get_parameters()->get_tau());
    }
    shift.resize(path->get_parameters()->get_ndim());
}

void Center_of_Mass::attempt(){
    Move_Base::attempt();
    
    path->get_beads()->set_old_data();
    
    for(fVector::iterator it = shift.begin(); it != shift.end(); it++ )
        *it = path->get_util()->randgaussian(delta);

    
    for(int slice = 0; slice < path->get_parameters()->get_num_timeslices(); slice++){
        old_action += pa->get_action(slice,0);
    }
    
    path->get_beads()->shift_all(ptcl, shift);
    
    for(int slice = 0; slice < path->get_parameters()->get_num_timeslices(); slice++){
        new_action += pa->get_action(slice,0);
    }
    
    if(check_move())
        accept();
    else
        reject();
}

void Center_of_Mass::accept(){
    Move_Base::accept();
    iVector lc = path->get_last_changed();
    lc.push_back(ptcl);
    path->set_last_changed(lc);
}

/************************************************************
 
                    BISECTION MOVE CLASS
 
 ************************************************************/

Bisection::Bisection(boost::shared_ptr<Path> path) : Move_Base(path){
    multistep_dist = path->get_multistep_dist();
}

void Bisection::attempt(){
    Move_Base::attempt();
    
    start = path->get_util()->randint(path->get_parameters()->get_num_timeslices());
    
    path->get_beads()->set_old_data();
    
    for(int a = 0; a < multistep_dist+1; a++){
        int slice = (start + a)%path->get_parameters()->get_num_timeslices();
        old_action += pa->get_action(slice,0);
    }
    
    level_move(ptcl, start, multistep_dist);
    
    for(int a = 0; a < multistep_dist+1; a++ ){
        int slice = (start + a)%path->get_parameters()->get_num_timeslices();
        new_action += pa->get_action(slice,0);
    }
    
    if(check_move())
        accept();
    else
        reject();
}

void Bisection::level_move(int ptcl, int start, int m){
    if(m != 1 && m%2 == 0){
        int slice = (start + m/2);
        float tau1 = (m/2)*path->get_parameters()->get_tau();
        fVector move(0);
        ffVector bds = path->get_beads()->get_pair_same_path(ptcl, start, m);
        fVector aved = path->get_util()->avedist(bds,path->get_parameters()->get_box_size());
        float width = sqrt(path->get_parameters()->get_lambda()*tau1);
        for(int ndim = 0; ndim < path->get_parameters()->get_ndim(); ndim++){
            move.push_back(aved[ndim] + path->get_util()->randgaussian(width));
        }
        path->get_beads()->set_bead_data(ptcl, slice, move);
        level_move(ptcl, start, m/2);
        level_move(ptcl, slice, m/2);
    }
}

void Bisection::accept(){
    
    iVector chd_ptcl;
    chd_ptcl.push_back(ptcl);
    
    path->set_last_changed(chd_ptcl);
    path->set_last_start_end(start, start+multistep_dist);
    path->get_beads()->set_prev_perm();
    
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
}

void Perm_Bisection::attempt(){
    Move_Base::attempt();

    start = path->get_util()->randint(path->get_parameters()->get_num_timeslices());
    
    int ls = path->get_last_start_end()[0];
    int le = path->get_last_start_end()[1];
    
    iVector lc;
    
    if(le > path->get_parameters()->get_num_timeslices() && start < le% path->get_parameters()->get_num_timeslices())
        lc = path->get_beads()->get_rep_swap_order();
    else
        lc = path->get_last_changed();
    
    if((start+multistep_dist > ls &&start+multistep_dist <=le) ||
                ((start >=ls) && (start <=le))|| (le > path->get_parameters()->get_num_timeslices()
                                                && start < le % path->get_parameters()->get_num_timeslices()))
        ptable->recalc_perms(lc, start);
    
    path->get_beads()->set_old_data();
    
    for(int a = 0; a < multistep_dist+1; a++){
        int slice = (start + a)%path->get_parameters()->get_num_timeslices();
        old_action += pa->get_action(slice,0);
    }
    
    iVector identity(path->get_beads()->get_num_particles());
    iota(identity.begin(),identity.end(),0);
    iVector chosenPerm = identity;
    iVector origpart(0);
    permed_parts = std::vector<int>(0);
    
    if(path->get_parameters()->is_boson()){
        chosenPerm = ptable->pick_permutation(0, start);
        
        for(iVector::iterator it = identity.begin(); it != identity.end(); it++)
            if(*it != chosenPerm[*it]){
                origpart.push_back(*it);
                permed_parts.push_back(chosenPerm[*it]);
            }
    }
    
    path->get_beads()->set_permutation(origpart, permed_parts, start, multistep_dist);
    path->get_beads()->permute();
    
    if(permed_parts.size() == 0){
        level_move(ptcl, start, multistep_dist);
        permed_parts.push_back(ptcl);
    }
    else
        for(iVector::iterator ptcl = origpart.begin(); ptcl !=origpart.end(); ptcl++)
            level_move(*ptcl, start, multistep_dist);
    
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
    path->set_last_changed(permed_parts);
    path->set_last_start_end(start, start + multistep_dist);
    path->get_beads()->set_prev_perm();
    
    Move_Base::accept();

}
void Perm_Bisection::reject(){
    path->get_beads()->permute(true);
    Bisection::reject();
    path->get_beads()->reset_permute();
    
    float check_action = 0;
    for(int a = 0; a < multistep_dist+1; a++ ){
        int slice = (start + a)%path->get_parameters()->get_num_timeslices();
        check_action += pa->get_action(slice,0);
    }
    
}



