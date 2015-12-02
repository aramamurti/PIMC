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
    ptcl = (int) path->get_util()->randnormed(path->get_beads()->get_num_particles())%path->get_beads()->get_num_particles();
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
    
    for(dVector::iterator it = shift.begin(); it != shift.end(); it++ )
        *it = path->get_util()->randgaussian(delta);

    
    for(int slice = 0; slice < path->get_parameters()->get_num_timeslices(); slice++){
        old_action += pa->get_action(slice,0, false);
    }
    
    changed_particles = path->get_beads()->shift_all(ptcl, shift);
    
    for(int slice = 0; slice < path->get_parameters()->get_num_timeslices(); slice++){
        new_action += pa->get_action(slice,0, false);
    }
    
    if(check_move())
        accept();
    else
        reject();
}

void Center_of_Mass::accept(){
    Move_Base::accept();
    iVector lc = path->get_last_changed();
    for(iVector::iterator it = changed_particles.begin(); it != changed_particles.end(); it++){
        lc.push_back(*it);
    }
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
        old_action += pa->get_action(slice,0, false);
    }
    
    level_move(ptcl, start, multistep_dist);
    
    for(int a = 0; a < multistep_dist+1; a++ ){
        int slice = (start + a)%path->get_parameters()->get_num_timeslices();
        new_action += pa->get_action(slice,0, false);
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
        dVector move(0);
        ddVector bds = path->get_beads()->get_pair_same_path(ptcl, start, m);
        iVector chd_ptcls = path->get_beads()->get_pair_rows(ptcl, start, m);
        changed_particles.push_back(chd_ptcls[0]);
        changed_particles.push_back(chd_ptcls[1]);
        dVector aved = path->get_util()->avedist(bds,path->get_parameters()->get_box_size());
        double width = sqrt(path->get_parameters()->get_lambda()*tau1);
        for(int ndim = 0; ndim < path->get_parameters()->get_ndim(); ndim++){
            move.push_back(aved[ndim] + path->get_util()->randgaussian(width));
        }
        path->get_beads()->set_bead_data(ptcl, slice, move);
        level_move(ptcl, start, m/2);
        level_move(ptcl, slice, m/2);
    }
}

void Bisection::accept(){
    
    iVector lc = path->get_last_changed();
    auto last = std::unique(changed_particles.begin(), changed_particles.end());
    changed_particles.erase(last, changed_particles.end());

    for(iVector::iterator it = changed_particles.begin(); it != changed_particles.end(); it++){
        lc.push_back(*it);
    }
    path->set_last_changed(lc);
    
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
    
    for(int slice = 0; slice < path->get_parameters()->get_num_timeslices(); slice++)
        ptable->recalc_perms(path->get_last_changed(), slice);
    
    path->get_beads()->set_old_data();
    
    for(int a = 0; a < multistep_dist+1; a++){
        int slice = (start + a)%path->get_parameters()->get_num_timeslices();
        old_action += pa->get_action(slice,0, false);
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
        new_action += pa->get_action(slice,0, false);
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
    
    Move_Base::accept();

}
void Perm_Bisection::reject(){
    path->get_beads()->permute(true);
    Bisection::reject();
    path->get_beads()->reset_permute();
    path->set_last_changed(iVector(0));
}

/***************************************************************************************
 
                                INSERT  WORM MOVE
 
 ***************************************************************************************/

Insert::Insert(boost::shared_ptr<Path> path) : Move_Base(path){

}

void Insert::attempt(){
    Move_Base::attempt();
    
    int m = path->get_util()->randint(path->get_parameters()->get_mbar())+1;
    int start_slice = path->get_util()->randint(path->get_parameters()->get_num_timeslices());
    
    dVector pos(0);
    for(int i = 0; i < path->get_parameters()->get_ndim();i++){
        pos.push_back(path->get_util()->randnormed(path->get_parameters()->get_box_size()));
    }
    
    path->get_beads()->initialize_worm();
    path->get_beads()->new_worm(pos, start_slice);
    
    for(int i = 1; i < m; i++){
        dVector new_pos(0);
        double width = sqrt(path->get_parameters()->get_lambda()*path->get_parameters()->get_tau());
        for(int ndim = 0; ndim < path->get_parameters()->get_ndim(); ndim++){
            new_pos.push_back(pos[ndim] + path->get_util()->randgaussian(width));
        }
        path->get_beads()->worm_push_back(new_pos);
        pos = new_pos;
    }
    
    mu_shift = path->get_parameters()->get_mu()*m;
    
    new_action = 0;
    for(int a = 0; a < m; a++ ){
        int slice = (start_slice + a)%path->get_parameters()->get_num_timeslices();
        new_action += pa->get_action(slice,0,true);
    }
    
    if(check_move())
        accept();
    else
        reject();

    
}
bool Insert::check_move(){
    if(path->get_util()->randnormed(1)<path->get_parameters()->get_C0()*exp(-new_action + mu_shift*path->get_parameters()->get_tau()))
        return true;
    else
        return false;
}
void Insert::accept(){
    Move_Base::accept();
    path->set_worm(true);
}

void Insert::reject(){
    path->get_beads()->remove_worm();
}


/***************************************************************************************
 
                                REMOVE WORM MOVE
 
 ***************************************************************************************/

Remove::Remove(boost::shared_ptr<Path> path) : Move_Base(path){
    
}

void Remove::attempt(){
    Move_Base::attempt();
    
    int worm_size =path->get_beads()->get_worm_size();
    if(worm_size > path->get_parameters()->get_mbar()){
        reject();
        return;
    }
    
    mu_shift = path->get_parameters()->get_mu()*worm_size;
    
    old_action = 0;
    for(int a = 0; a < worm_size; a++ ){
        int slice = (path->get_beads()->get_worm_indices()[0].second + a)%path->get_parameters()->get_num_timeslices();
        old_action += pa->get_action(slice,0,true);
    }
    
    if(check_move())
        accept();
    else
        reject();

    
}
bool Remove::check_move(){
    if(path->get_util()->randnormed(1)<path->get_parameters()->get_C0()*exp((old_action - mu_shift*path->get_parameters()->get_tau())))
        return true;
    else
        return false;
}
void Remove::accept(){
    Move_Base::accept();
    path->get_beads()->remove_worm();
    path->set_worm(false);
}

void Remove::reject(){}



/***************************************************************************************
 
                                ADVANCE HEAD MOVE
 
 ***************************************************************************************/


Advance_Head::Advance_Head(boost::shared_ptr<Path> path) : Move_Base(path){
    
}

void Advance_Head::attempt(){
    Move_Base::attempt();
    
    int m = path->get_util()->randint(path->get_parameters()->get_mbar())+1;
    int start_slice = path->get_beads()->get_worm_indices()[0].second;
    
    dVector pos = path->get_beads()->get_worm_bead_data(0, start_slice);
    for(int i = 0; i < m; i++){
        dVector new_pos(0);
        double width = sqrt(path->get_parameters()->get_lambda()*path->get_parameters()->get_tau());
        for(int ndim = 0; ndim < path->get_parameters()->get_ndim(); ndim++){
            new_pos.push_back(pos[ndim] + path->get_util()->randgaussian(width));
        }
        path->get_beads()->worm_push_front(new_pos);
        pos = new_pos;
    }
    
    mu_shift = path->get_parameters()->get_mu()*m;
    
    new_action = 0;
    
    for(int a = 1; a <= m; a++){
        int slice = (start_slice - a+path->get_parameters()->get_num_timeslices())%path->get_parameters()->get_num_timeslices();
        new_action += pa->get_action(slice,0,true);
    }
    
    if(check_move())
        accept();
    else
        reject(m);
    
    
}
bool Advance_Head::check_move(){
    if(path->get_util()->randnormed(1)<exp(-new_action + mu_shift*path->get_parameters()->get_tau()))
        return true;
    else
        return false;
}
void Advance_Head::accept(){
    Move_Base::accept();
}

void Advance_Head::reject(int m){
    for(int i = 0; i < m; i++){
        path->get_beads()->worm_pop_front(true);
    }
}



/***************************************************************************************
 
                                RECEDE HEAD MOVE
 
 ***************************************************************************************/

Recede_Head::Recede_Head(boost::shared_ptr<Path> path) : Move_Base(path){
    
}

void Recede_Head::attempt(){
    Move_Base::attempt();
    
    int m = path->get_util()->randint(path->get_parameters()->get_mbar())+1;
    int start_slice = path->get_beads()->get_worm_indices()[0].second;

    
    int worm_size =path->get_beads()->get_worm_size();
    if(worm_size < m){
        reject();
        return;
    }
    
    mu_shift = path->get_parameters()->get_mu()*worm_size;
    
    old_action = 0;
    for(int a = 0; a < m; a++ ){
        int slice = (start_slice + a)%path->get_parameters()->get_num_timeslices();
        old_action += pa->get_action(slice,0,true,m);
    }
    
    if(check_move())
        accept(m);
    else
        reject();

    
    
}
bool Recede_Head::check_move(){
    if(path->get_util()->randnormed(1)<exp(old_action - mu_shift*path->get_parameters()->get_tau()))
        return true;
    else
        return false;
}
void Recede_Head::accept(int m){
    Move_Base::accept();
    for(int i = 0; i < m; i++){
        path->get_beads()->worm_pop_front(true);
    }
    if(path->get_beads()->get_worm_size() == 0)
        path->set_worm(false);
}

void Recede_Head::reject(){

}




/***************************************************************************************
 
                                ADVANCE TAIL MOVE
 
 ***************************************************************************************/

Advance_Tail::Advance_Tail(boost::shared_ptr<Path> path) : Move_Base(path){
    
}

void Advance_Tail::attempt(){
    Move_Base::attempt();
    
    int m = path->get_util()->randint(path->get_parameters()->get_mbar())+1;
    int start_slice = path->get_beads()->get_worm_indices()[1].second;
    
    dVector pos = path->get_beads()->get_worm_bead_data(path->get_beads()->get_worm_indices()[1].first, start_slice);
    for(int i = 0; i < m; i++){
        dVector new_pos(0);
        double width = sqrt(path->get_parameters()->get_lambda()*path->get_parameters()->get_tau());
        for(int ndim = 0; ndim < path->get_parameters()->get_ndim(); ndim++){
            new_pos.push_back(pos[ndim] + path->get_util()->randgaussian(width));
        }
        path->get_beads()->worm_push_back(new_pos);
        pos = new_pos;
    }
    
    mu_shift = path->get_parameters()->get_mu()*m;
    
    new_action = 0;
    
    for(int a = 1; a <= m; a++){
        int slice = (start_slice + a)%path->get_parameters()->get_num_timeslices();
        new_action += pa->get_action(slice,0,true);
    }
    
    if(check_move())
        accept();
    else
        reject(m);
    
    
}
bool Advance_Tail::check_move(){
    if(path->get_util()->randnormed(1)<exp(-new_action + mu_shift*path->get_parameters()->get_tau()))
        return true;
    else
        return false;
}
void Advance_Tail::accept(){
    Move_Base::accept();
}

void Advance_Tail::reject(int m){
    for(int i = 0; i < m; i++){
        path->get_beads()->worm_pop_front(true);
    }
}



/***************************************************************************************
 
                                RECEDE TAIL MOVE
 
 ***************************************************************************************/
Recede_Tail::Recede_Tail(boost::shared_ptr<Path> path) : Move_Base(path){
    
}

void Recede_Tail::attempt(){
    Move_Base::attempt();
    
    int m = path->get_util()->randint(path->get_parameters()->get_mbar())+1;
    int start_slice = path->get_beads()->get_worm_indices()[1].second;
    
    int worm_size =path->get_beads()->get_worm_size();
    if(worm_size < m){
        reject();
        return;
    }
    
    mu_shift = path->get_parameters()->get_mu()*worm_size;
    
    old_action = 0;
    for(int a = 0; a < m; a++ ){
        int slice = (start_slice - a + path->get_parameters()->get_num_timeslices())%path->get_parameters()->get_num_timeslices();
        old_action += pa->get_action(slice,0,true,0,m);
    }
    
    if(check_move())
        accept(m);
    else
        reject();
    
}
bool Recede_Tail::check_move(){
    if(path->get_util()->randnormed(1)<exp(old_action - mu_shift*path->get_parameters()->get_tau()))
        return true;
    else
        return false;
}
void Recede_Tail::accept(int m){
    Move_Base::accept();
    for(int i = 0; i < m; i++){
        path->get_beads()->worm_pop_back(true);
    }
    if(path->get_beads()->get_worm_size() == 0)
        path->set_worm(false);
}

void Recede_Tail::reject(){
    
}

/*******************************************************************
 
                            OPEN MOVE
 
 ******************************************************************/



/*******************************************************************
 
                            CLOSE MOVE
 
 ******************************************************************/


/*******************************************************************
 
                            SWAP HEAD MOVE
 
 ******************************************************************/


/*******************************************************************
 
                            SWAP TAIL MOVE
 
 ******************************************************************/
