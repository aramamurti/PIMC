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
    ka = boost::shared_ptr<Kinetic_Action>(new Kinetic_Action(this->path));
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
    worm_move = false;
    worm_nec = false;
    
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
    
    path->get_beads()->shift_all(ptcl, shift);
    
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
        iVector lc = path->get_last_changed();
        for(iVector::iterator it = changed_particles.begin(); it != changed_particles.end(); it++){
            lc.push_back(*it);
        }
        auto last = std::unique(lc.begin(), lc.end());
        lc.erase(last, lc.end());
        path->set_last_changed(lc);
    }
    Move_Base::accept();
    path->put_in_box(changed_particles,0,path->get_parameters()->get_num_timeslices());

}

/************************************************************
 
 BISECTION MOVE CLASS
 
 ************************************************************/

Bisection::Bisection(boost::shared_ptr<Path> path) : Move_Base(path){
    multistep_dist = path->get_multistep_dist();
    worm_move = false;
    worm_nec = false;
    
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
    std::sort(changed_particles.begin(), changed_particles.end());
    auto last = std::unique(changed_particles.begin(), changed_particles.end());
    changed_particles.erase(last, changed_particles.end());
    
    if(path->get_parameters()->get_move_list()[2]){
        iVector lc = path->get_last_changed();
        
        for(iVector::iterator it = changed_particles.begin(); it != changed_particles.end(); it++){
            lc.push_back(*it);
        }
        last = std::unique(lc.begin(), lc.end());
        lc.erase(last, lc.end());
        path->set_last_changed(lc);
        
        path->set_last_start_end(start, start+multistep_dist);
        path->get_beads()->set_prev_perm();
    }
    Move_Base::accept();
    path->put_in_box(changed_particles, start, start+multistep_dist);
}

void Bisection::reject(){
    Move_Base::reject();
}

/**********************************************************************************
 
 PERMUTATION BISECTION MOVE CLASS
 
 ***********************************************************************************/

Perm_Bisection::Perm_Bisection(boost::shared_ptr<Path> path) : Bisection(path){
    ptable = boost::shared_ptr<Permutation_Table>(new Permutation_Table(path));
    worm_move = false;
    worm_nec = false;
    
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
        old_action += pa->get_action(slice,0, false);
    }
    
    iVector identity(num_particles);
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
    path->get_beads()->reset_permute();
    path->set_last_changed(iVector(0));
    Bisection::reject();
}

/***************************************************************************************
 
 INSERT  WORM MOVE
 
 ***************************************************************************************/

Insert::Insert(boost::shared_ptr<Path> path) : Move_Base(path){
    
    worm_move = true;
    worm_nec = false;
    
    move_name = "Insert";
    
}

void Insert::attempt(){
    Move_Base::attempt();
    
    // std::cout << "Insert Attempt" << std::endl;
    
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
 //       path->get_beads()->check_key_map();
        pos = new_pos;
    }
    
    mu_shift = path->get_parameters()->get_mu()*m;
    
    path->set_worm(true);
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
    path->set_worm(true);
    Move_Base::accept();
    path->put_worm_in_box();


}

void Insert::reject(){
     // std::cout << "Rejected" << std::endl;
    path->get_beads()->remove_worm();
    path->set_worm(false);
    Move_Base::reject();
    
}


/***************************************************************************************
 
 REMOVE WORM MOVE
 
 ***************************************************************************************/

Remove::Remove(boost::shared_ptr<Path> path) : Move_Base(path){
    worm_move = true;
    worm_nec = true;
    
    move_name = "Remove";
    
}

void Remove::attempt(){
    Move_Base::attempt();
    
    //// std::cout << "Remove Attempt" << std::endl;
    
    
    int worm_size =path->get_beads()->get_worm_size();
    if(worm_size > path->get_parameters()->get_mbar()){
        reject();
        return;
    }
    
    mu_shift = path->get_parameters()->get_mu()*worm_size;
    
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
    if(path->get_util()->randnormed(1)<path->get_parameters()->get_C0()*exp(old_action - mu_shift*path->get_parameters()->get_tau()))
        return true;
    else
        return false;
}
void Remove::accept(){
     // std::cout << "Accepted" << std::endl;
    
    path->get_beads()->remove_worm();
    path->set_worm(false);
    Move_Base::accept();

}

void Remove::reject(){
    Move_Base::reject();
}


/***************************************************************************************
 
 ADVANCE HEAD MOVE
 
 ***************************************************************************************/


Advance_Head::Advance_Head(boost::shared_ptr<Path> path) : Move_Base(path){
    worm_move = true;
    worm_nec = true;
    
    move_name = "Advance Head";
    
}

void Advance_Head::attempt(){
    Move_Base::attempt();
    
    start_end.resize(0);
    // std::cout << "Advance Head Attempt" << std::endl;
    
    
    int m = path->get_util()->randint(path->get_parameters()->get_mbar())+1;
    
    std::pair<int, int> head_index = path->get_beads()->get_worm_indices()[0];
    int start_slice = head_index.second;
    
    for(int a = 1; a <= m; a++){
        int slice = (start_slice - a+path->get_parameters()->get_num_timeslices())%path->get_parameters()->get_num_timeslices();
        old_action += pa->get_action(slice,0,true);
    }
    
    std::vector<std::pair<int, int> > ht = path->get_beads()->get_worm_indices();
    
    //    int dist;
    //    if(ht[1].second >= ht[0].second){
    //        dist = path->get_parameters()->get_num_timeslices()-(ht[1].second-ht[0].second) -1;
    //    }
    //    else{
    //        dist = ht[1].second-ht[0].second-1;
    //    }
    //
    //    if(dist <= m){
    //        dist += path->get_parameters()->get_num_timeslices();
    //    }
    
    dVector pos = path->get_beads()->get_worm_bead_data(0, start_slice);
    dVector end_pos = path->get_beads()->get_worm_bead_data(ht[1].first, ht[1].second);
    
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
    
    
    for(int a = 1; a <= m; a++){
        int slice = (start_slice - a+path->get_parameters()->get_num_timeslices())%path->get_parameters()->get_num_timeslices();
        new_action += pa->get_action(slice,0,true);
    }
    
    if(check_move()){
        std::pair<int, int> new_head_index = path->get_beads()->get_worm_indices()[0];
        if(head_index.second > new_head_index.second){
            changed_particles.push_back(head_index.first);
            start_end.push_back(std::pair<int,int>(new_head_index.second, head_index.second));
        }
        else{
            changed_particles.push_back(0);
            changed_particles.push_back(1);
            start_end.push_back(std::pair<int,int>(new_head_index.second,path->get_parameters()->get_num_timeslices()-1));
            start_end.push_back(std::pair<int,int>(0, head_index.second));
            
        }
        accept();
    }
    else
        reject(m);
    
    
}
bool Advance_Head::check_move(){
    if(path->get_util()->randnormed(1)<exp(-(new_action-old_action) + mu_shift*path->get_parameters()->get_tau()))
        return true;
    else
        return false;
}
void Advance_Head::accept(){
    Move_Base::accept();
    path->put_worm_in_box(changed_particles,start_end);
}

void Advance_Head::reject(int m){
     // std::cout << "Rejected" << std::endl;
    
    for(int i = 0; i < m; i++){
        path->get_beads()->worm_pop_front(true);
    }
    Move_Base::reject();
}



/***************************************************************************************
 
 RECEDE HEAD MOVE
 
 ***************************************************************************************/

Recede_Head::Recede_Head(boost::shared_ptr<Path> path) : Move_Base(path){
    worm_move = true;
    worm_nec = true;
    
    move_name = "Recede Head";
}

void Recede_Head::attempt(){
    Move_Base::attempt();
    
    //// std::cout << "Recede Head Attempt" << std::endl;
    
    
    int m = path->get_util()->randint(path->get_parameters()->get_mbar())+1;
    int start_slice = path->get_beads()->get_worm_indices()[0].second;
    
    
    int worm_size =path->get_beads()->get_worm_size();
    if(worm_size < m){
        reject();
        return;
    }
    
    mu_shift = path->get_parameters()->get_mu()*m;
    
    for(int a = 0; a < m; a++ ){
        int slice = (start_slice + a)%path->get_parameters()->get_num_timeslices();
        old_action += pa->get_action(slice,0,true);
    }
    
    for(int a = 0; a < m; a++ ){
        int slice = (start_slice + a)%path->get_parameters()->get_num_timeslices();
        new_action += pa->get_action(slice,0,true,m);
    }
    
    
    if(check_move())
        accept(m);
    else
        reject();
    
    
    
}
bool Recede_Head::check_move(){
    if(path->get_util()->randnormed(1)<exp(-(new_action-old_action) - mu_shift*path->get_parameters()->get_tau()))
        return true;
    else
        return false;
}
void Recede_Head::accept(int m){
     // std::cout << "Accepted" << std::endl;
    
    for(int i = 0; i < m; i++){
        path->get_beads()->worm_pop_front(true);
    }
    
    if(path->get_beads()->get_worm_size() == 0)
        path->set_worm(false);

    Move_Base::accept();

}

void Recede_Head::reject(){
    Move_Base::reject();
}



/***************************************************************************************
 
 ADVANCE TAIL MOVE
 
 ***************************************************************************************/

Advance_Tail::Advance_Tail(boost::shared_ptr<Path> path) : Move_Base(path){
    worm_move = true;
    worm_nec = true;
    
    move_name = "Advance Tail";
}

void Advance_Tail::attempt(){
    Move_Base::attempt();
    start_end.resize(0);
    
    // std::cout << "Advance Tail Attempt" << std::endl;
    
    
    int m = path->get_util()->randint(path->get_parameters()->get_mbar())+1;
    std::pair<int,int> tail_index = path->get_beads()->get_worm_indices()[1];
    int start_slice = tail_index.second;
    
    for(int a = 1; a <= m; a++){
        int slice = (start_slice + a)%path->get_parameters()->get_num_timeslices();
        old_action += pa->get_action(slice,0,true);
    }
    
    std::vector<std::pair<int, int> > ht = path->get_beads()->get_worm_indices();
    
    
    
    dVector pos = path->get_beads()->get_worm_bead_data(path->get_beads()->get_worm_indices()[1].first, start_slice);
    dVector end_pos = path->get_beads()->get_worm_bead_data(ht[0].first, ht[0].second);
    
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
    
    for(int a = 1; a <= m; a++){
        int slice = (start_slice + a)%path->get_parameters()->get_num_timeslices();
        new_action += pa->get_action(slice,0,true);
    }
    
    if(check_move()){
        std::pair<int,int> new_tail_index = path->get_beads()->get_worm_indices()[1];
        if(new_tail_index.first == tail_index.first){
            changed_particles.push_back(tail_index.first);
            start_end.push_back(std::pair<int,int>(tail_index.second,new_tail_index.second));
        }
        else{
            changed_particles.push_back(tail_index.first);
            changed_particles.push_back(new_tail_index.first);
            start_end.push_back(std::pair<int,int>(tail_index.second,path->get_parameters()->get_num_timeslices()-1));
            start_end.push_back(std::pair<int,int>(0,new_tail_index.second));
        }
        accept();
    }
    else
        reject(m);
    
    
}
bool Advance_Tail::check_move(){
    if(path->get_util()->randnormed(1)<exp(-(new_action-old_action) + mu_shift*path->get_parameters()->get_tau()))
        return true;
    else
        return false;
}
void Advance_Tail::accept(){
    Move_Base::accept();
    path->put_worm_in_box(changed_particles, start_end);
}

void Advance_Tail::reject(int m){
     // std::cout << "Rejected" << std::endl;
    
    for(int i = 0; i < m; i++){
        path->get_beads()->worm_pop_back(true);
    }
    Move_Base::reject();

}



/***************************************************************************************
 
 RECEDE TAIL MOVE
 
 ***************************************************************************************/
Recede_Tail::Recede_Tail(boost::shared_ptr<Path> path) : Move_Base(path){
    worm_move = true;
    worm_nec = true;
    
    move_name = "Recede Tail";
}

void Recede_Tail::attempt(){
    Move_Base::attempt();
    
     // std::cout << "Recede Tail Attempt" << std::endl;
    
    int m = path->get_util()->randint(path->get_parameters()->get_mbar())+1;
    int start_slice = path->get_beads()->get_worm_indices()[1].second;
    
    int worm_size =path->get_beads()->get_worm_size();
    if(worm_size < m){
        reject();
        return;
    }
    
    mu_shift = path->get_parameters()->get_mu()*m;
    
    for(int a = 0; a < m; a++ ){
        int slice = (start_slice - a + path->get_parameters()->get_num_timeslices())%path->get_parameters()->get_num_timeslices();
        old_action += pa->get_action(slice,0,true);
    }
    
    for(int a = 0; a < m; a++ ){
        int slice = (start_slice - a + path->get_parameters()->get_num_timeslices())%path->get_parameters()->get_num_timeslices();
        new_action += pa->get_action(slice,0,true,0,m);
    }
    
    if(check_move())
        accept(m);
    else
        reject();
    
}
bool Recede_Tail::check_move(){
    if(path->get_util()->randnormed(1)<exp(-(new_action-old_action) - mu_shift*path->get_parameters()->get_tau()))
        return true;
    else
        return false;
}
void Recede_Tail::accept(int m){
     // std::cout << "Accepted" << std::endl;
    
    for(int i = 0; i < m; i++){
        path->get_beads()->worm_pop_back(true);
    }
    if(path->get_beads()->get_worm_size() == 0)
        path->set_worm(false);

    Move_Base::accept();
}

void Recede_Tail::reject(){
    Move_Base::reject();
}

/*******************************************************************
 
 OPEN MOVE
 
 ******************************************************************/

Open::Open(boost::shared_ptr<Path> path) : Move_Base(path){
    worm_move = true;
    worm_nec = false;
    
    move_name = "Open";
    
}

void Open::attempt(){
    Move_Base::attempt();
    
    if(num_particles == 0){
        reject();
        return;
    }
    
    // std::cout << "Open Attempt" << std::endl;
    
    m = path->get_util()->randint(path->get_parameters()->get_mbar())+1;
    start_slice = path->get_util()->randint(path->get_parameters()->get_num_timeslices());
    
    mu_shift = path->get_parameters()->get_mu()*m;
    changed_particles = path->get_beads()->get_changed_ptcls(ptcl);
    
    for(int i = 1; i < m; i ++){
        int slice = (start_slice+i);
        int row = ptcl;
        if(slice >= path->get_parameters()->get_num_timeslices()){
            slice = slice%path->get_parameters()->get_num_timeslices();
            row = path->get_beads()->get_pair_rows(ptcl, start_slice, m)[1];
        }
        old_action += pa->get_action_single_particle(row, slice);
    }
    
    if(check_move())
        accept();
    else
        reject();
    
}

bool Open::check_move(){
    double norm = path->get_parameters()->get_C0()*(num_particles-changed_particles.size())/pow(path->get_parameters()->get_box_size(),3);
    double rho0 = 1.0/pow((4*M_PI*path->get_parameters()->get_lambda()*m*path->get_parameters()->get_tau()),path->get_parameters()->get_ndim()/2.)*exp(-ka->get_action_single_particle(ptcl, start_slice, m));
    if(path->get_util()->randnormed(1)<(norm/rho0*exp(old_action - mu_shift*path->get_parameters()->get_tau())))
        return true;
    else
        return false;
}

void Open::accept(){
     // std::cout << "Accepted" << std::endl;
    int row = path->get_beads()->get_pair_rows(ptcl, start_slice, m)[1];
    rows_to_remove = path->get_beads()->move_path_to_worm(row, (start_slice+m)%path->get_parameters()->get_num_timeslices(), m-1);
    path->get_beads()->remove_old_list_references(rows_to_remove);
    path->set_worm(true);
    Move_Base::accept();
}

/*******************************************************************
 
 CLOSE MOVE
 
 ******************************************************************/

Close::Close(boost::shared_ptr<Path> path) : Move_Base(path){
    worm_move = true;
    worm_nec = true;
    
    move_name = "Close";
}

void Close::attempt(){
    Move_Base::attempt();
    start_end.resize(0);
    
    
    // std::cout << "Close Attempt" << std::endl;
    
    ht = path->get_beads()->get_worm_indices();
    
    int old_head_col = ht[0].second;
    int old_tail_col = ht[1].second;
    if(old_tail_col >= old_head_col){
        m = path->get_parameters()->get_num_timeslices()-(old_tail_col-old_head_col) -1;
    }
    else{
        m = old_head_col-old_tail_col-1;
    }
    
    if(m!=0)
        for(int slice = old_tail_col+1; ; slice++){
            slice = slice%path->get_parameters()->get_num_timeslices();
            if(slice == old_head_col)
                break;
            old_action += pa->get_action(slice, 0, true);
        }
    
    mu_shift = m*path->get_parameters()->get_mu();
    
    if(m > path->get_parameters()->get_mbar()){
        reject();
        return;
    }
    path->get_beads()->temporarily_close_worm();
    if(m != 0){
        
        
        ddVector new_pos(m, dVector(path->get_parameters()->get_ndim(),0));
        dVector start_pos = path->get_beads()->get_worm_bead_data(ht[1].first, ht[1].second);
        dVector end_pos = path->get_beads()->get_worm_bead_data(ht[0].first, ht[0].second);
        for(int n = 0; n < m; n++){
            double tau1 = (m-n)*path->get_parameters()->get_tau();
            dVector avex(path->get_parameters()->get_ndim());
            for(dVector::iterator it = avex.begin(); it!= avex.end(); it++){
                *it = (tau1*start_pos[it-avex.begin()]+path->get_parameters()->get_tau()*end_pos[it-avex.begin()])/(tau1+path->get_parameters()->get_tau());
                double box_size = path->get_parameters()->get_box_size();
                if(box_size != -1)
                    *it = path->get_util()->per_bound_cond(*it, box_size);
            }
            double width = sqrt(2 * path->get_parameters()->get_lambda()/(1/path->get_parameters()->get_tau()+tau1));
            dVector bead_pos(path->get_parameters()->get_ndim());
            for(dVector::iterator it = bead_pos.begin(); it != bead_pos.end(); it++){
                *it = avex[it-bead_pos.begin()] + path->get_util()->randgaussian(width);
            }
            new_pos[n]=bead_pos;
            start_pos = bead_pos;
        }
        
        if(old_tail_col >= old_head_col){
            int counter = 0;
            int worm_row = ht[1].first;
            int worm_col = old_tail_col+1;
            changed_particles.push_back(worm_row);
            start_end.push_back(std::pair<int, int>(worm_col, path->get_parameters()->get_num_timeslices()-1));
            while(counter < m){
                if(worm_col >= path->get_beads()->get_worm_dims()[1]){
                    worm_row = 0;
                    worm_col = 0;
                    changed_particles.push_back(worm_row);
                }
                path->get_beads()->set_worm_bead_data(worm_row, worm_col, new_pos[counter]);
                worm_col++;
                counter++;
            }
            start_end.push_back(std::pair<int,int>(0, worm_col-1));
        }
        else{
            int counter = 0;
            int worm_row = 0;
            int worm_col = old_head_col-1;
            changed_particles.push_back(0);
            while(counter < m){
                path->get_beads()->set_worm_bead_data(worm_row, worm_col, new_pos[m-1-counter]);
                counter++;
                worm_col--;
            }
            start_end.push_back(std::pair<int, int>(worm_col+1,old_head_col-1));
            
        }
        
        for(int slice = old_tail_col+1; slice != old_head_col; slice++){
            slice = slice%path->get_parameters()->get_num_timeslices();
            if(slice == old_head_col)
                break;
            new_action += pa->get_action(slice, 0, true);
        }
    }
    
    if(check_move())
        accept();
    else
        reject();
    
}

bool Close::check_move(){
    double norm = path->get_parameters()->get_C0()*(num_particles+path->get_beads()->get_worm_indices()[1].second+1)/pow(path->get_parameters()->get_box_size(),3);
    double rho0 = 1.0/pow((4*M_PI*path->get_parameters()->get_lambda()*m*path->get_parameters()->get_tau()),path->get_parameters()->get_ndim()/2.)*exp(-ka->get_action_worm_head_tail(ht[0].second, ht[1].second, m));
    if(path->get_util()->randnormed(1)<(rho0/norm*exp(-(new_action-old_action) + mu_shift*path->get_parameters()->get_tau())))
        return true;
    else
        return false;
}

void Close::reject(){
     // std::cout << "Rejected" << std::endl;
    path->get_beads()->reopen_temp_closed_worm();
    Move_Base::reject();

}

void Close::accept(){
    
    path->put_worm_in_box(changed_particles, start_end);
    path->get_beads()->move_worm_to_path();
    path->get_beads()->remove_old_worm_references();

    path->set_worm(false);
    Move_Base::accept();


    
}

/*******************************************************************
 
 SWAP HEAD MOVE
 
 ******************************************************************/

Swap_Head::Swap_Head(boost::shared_ptr<Path> path) : Move_Base(path){
    worm_move = true;
    worm_nec = true;
    
    move_name = "Swap Head";
}

void Swap_Head::attempt(){
    
    Move_Base::attempt();
    start_end.resize(0);
    
    if(num_particles == 0){
        m = 0;
        reject();
        return;
    }
    
    
    // std::cout << "Swap Head Attempt" << std::endl;
    std::vector<std::pair<int, int> > ht = path->get_beads()->get_worm_indices();
    int head_row = ht[0].first;
    int head_col = ht[0].second;
    
    m = path->get_parameters()->get_mbar();
    
    
    std::vector<std::pair<size_t, dVector> > L_I = path->get_beads()->get_worm_head_neighbors(m);
    std::vector<std::pair<size_t, double> > rho0s_I;
    sig_I = 0;
    for(std::vector<std::pair<size_t, dVector> >::iterator it = L_I.begin(); it != L_I.end(); it++){
        ddVector pair;
        pair.push_back(path->get_beads()->get_worm_bead_data(head_row, head_col));
        pair.push_back(it->second);
        double rho0 = 1.0/pow((4*M_PI*path->get_parameters()->get_lambda()*path->get_parameters()->get_mbar()*path->get_parameters()->get_tau()),path->get_parameters()->get_ndim()/2.)*exp(-ka->get_action_pos(pair, path->get_parameters()->get_mbar()));
        rho0s_I.push_back(std::pair<size_t, double>(it->first, rho0));
        sig_I += rho0;
    }
    if(isnan(sig_I)){
        sig_I = 0;
        for(std::vector<std::pair<size_t, dVector> >::iterator it = L_I.begin(); it != L_I.end(); it++){
            ddVector pair;
            pair.push_back(path->get_beads()->get_worm_bead_data(head_row, head_col));
            pair.push_back(it->second);
            double rho0 = 1.0/pow((4*M_PI*path->get_parameters()->get_lambda()*path->get_parameters()->get_mbar()*path->get_parameters()->get_tau()),path->get_parameters()->get_ndim()/2.)*exp(-ka->get_action_pos(pair, path->get_parameters()->get_mbar()));
            rho0s_I.push_back(std::pair<size_t, double>(it->first, rho0));
            sig_I += rho0;
        }
        
    }
    
    double rn = path->get_util()->randnormed(sig_I);
    size_t choice = 0;
    double accumulate = 0;
    for(std::vector<std::pair<size_t, double> >::iterator it = rho0s_I.begin(); it != rho0s_I.end(); it++){
        accumulate += it->second;
        if(rn < accumulate){
            choice = it->first;
            break;
        }
    }
    
    if(choice == 0){
        m = 0;
        reject();
        return;
    }
    
    ptcl = path->get_beads()->get_bead_indices(choice).first;
    changed_particles.push_back(ptcl);
    
    if(head_col - path->get_parameters()->get_mbar() < 0){
        changed_particles = path->get_beads()->get_pair_rows(ptcl, head_col, path->get_parameters()->get_mbar());
    }
    std::sort(changed_particles.begin(), changed_particles.end());
    auto last = std::unique(changed_particles.begin(), changed_particles.end());
    changed_particles.erase(last, changed_particles.end());
    
    
    for(int a = 0; a < path->get_parameters()->get_mbar(); a++){
        int slice = (head_col - a + path->get_parameters()->get_num_timeslices())%path->get_parameters()->get_num_timeslices();
        if(changed_particles.size() == 1)
            old_action += pa->get_action_single_particle(ptcl, slice, true);
        else
            old_action += pa->get_action_multiple_particles(changed_particles, slice, true);
        old_action += pa->get_action(slice,0,true);
    }
    changed_particles.resize(0);
    
    
    size_t ksi_key = path->get_beads()->get_next_bead_key(choice, m);
    if(!path->get_beads()->check_if_neighbor(path->get_beads()->get_worm_key(head_row, head_col), ksi_key)){
        m = 0;
        reject();
        return;
    }
    
    
    std::vector<std::pair<size_t, dVector> > L_ksi = path->get_beads()->get_neighbors(ksi_key, -m);
    
    sig_ksi = 0;
    for(std::vector<std::pair<size_t, dVector> >::iterator it = L_I.begin(); it != L_I.end(); it++){
        ddVector pair;
        pair.push_back(path->get_beads()->get_bead_data(path->get_beads()->get_bead_indices(ksi_key).first, path->get_beads()->get_bead_indices(ksi_key).second));
        pair.push_back(it->second);
        double rho0 = 1.0/pow((4*M_PI*path->get_parameters()->get_lambda()*path->get_parameters()->get_mbar()*path->get_parameters()->get_tau()),path->get_parameters()->get_ndim()/2.)*exp(-ka->get_action_pos(pair, m));
        sig_ksi += rho0;
    }
    if(isnan(sig_ksi)){
        sig_ksi = 0;
        for(std::vector<std::pair<size_t, dVector> >::iterator it = L_I.begin(); it != L_I.end(); it++){
            ddVector pair;
            pair.push_back(path->get_beads()->get_bead_data(path->get_beads()->get_bead_indices(ksi_key).first, path->get_beads()->get_bead_indices(ksi_key).second));
            pair.push_back(it->second);
            double rho0 = 1.0/pow((4*M_PI*path->get_parameters()->get_lambda()*path->get_parameters()->get_mbar()*path->get_parameters()->get_tau()),path->get_parameters()->get_ndim()/2.)*exp(-ka->get_action_pos(pair, m));
            sig_ksi += rho0;
        }
    }
    
    
    rows_to_remove = path->get_beads()->swap_into_worm_head(path->get_beads()->get_bead_indices(choice).first, m);
    
    
    ddVector new_pos(m-1, dVector(path->get_parameters()->get_ndim(),0));
    
    
    int shift = path->get_beads()->get_worm_indices()[1].first - ht[1].first;
    int start_col = ht[0].second;
    int start_row = ht[0].first+shift;
    
    dVector start_pos = path->get_beads()->get_worm_bead_data(start_row, start_col);
    
    int end_col = start_col - m;
    int end_row = start_row;
    
    if(end_col < 0){
        end_col = end_col + path->get_parameters()->get_num_timeslices();
        end_row--;
    }
    
    if(start_row == end_row){
        changed_particles.push_back(start_row);
        start_end.push_back(std::pair<int,int>(end_col+1, start_col-1));
    }
    else{
        if(start_col != 0){
            changed_particles.push_back(start_row);
            start_end.push_back(std::pair<int,int>(0,start_col-1));
        }
        if(end_col != path->get_parameters()->get_num_timeslices()-1){
            changed_particles.push_back(end_row);
            start_end.push_back(std::pair<int,int>(end_col+1,path->get_parameters()->get_num_timeslices()-1));
        }
        
    }
    dVector end_pos = path->get_beads()->get_worm_bead_data(end_row, end_col);
    for(int n = 0; n < m-1; n++){
        double tau1 = (m-1-n)*path->get_parameters()->get_tau();
        dVector avex(path->get_parameters()->get_ndim());
        for(dVector::iterator it = avex.begin(); it!= avex.end(); it++){
            *it = (tau1*start_pos[it-avex.begin()]+path->get_parameters()->get_tau()*end_pos[it-avex.begin()])/(tau1+path->get_parameters()->get_tau());
            double box_size = path->get_parameters()->get_box_size();
            if(box_size != -1)
                *it = path->get_util()->per_bound_cond(*it, box_size);
        }
        double width = sqrt(2. * path->get_parameters()->get_lambda()/(1.0/path->get_parameters()->get_tau()+tau1));
        dVector bead_pos(path->get_parameters()->get_ndim());
        for(dVector::iterator it = bead_pos.begin(); it != bead_pos.end(); it++){
            *it = avex[it-bead_pos.begin()] + path->get_util()->randgaussian(width);
        }
        new_pos[n] = bead_pos;
        start_pos = bead_pos;
    }
    
    int counter = 0;
    int moving_head_col = start_col-1;
    int moving_head_row = start_row;
    while(counter < m-1){
        if(moving_head_col < 0){
            moving_head_row--;
            moving_head_col = path->get_parameters()->get_num_timeslices()-1;
        }
        path->get_beads()->set_worm_bead_data(moving_head_row, moving_head_col, new_pos[m-2-counter]);
        moving_head_col--;
        counter++;
    }
    
    for(int a = 0; a < path->get_parameters()->get_mbar(); a++){
        int slice = (start_col - a+path->get_parameters()->get_num_timeslices())%path->get_parameters()->get_num_timeslices();
        new_action += pa->get_action(slice,0,true);
    }
    
    
    if(isnan(new_action)){
        new_action = 0;
        for(int a = 0; a < path->get_parameters()->get_mbar(); a++){
            int slice = (start_col - a+path->get_parameters()->get_num_timeslices())%path->get_parameters()->get_num_timeslices();
            new_action += pa->get_action(slice,0, false);
        }
    }
    
    if(check_move())
        accept();
    else
        reject();
    
}

bool Swap_Head::check_move(){
    if(path->get_util()->randnormed(1)< sig_I/sig_ksi * exp(-(new_action-old_action)))
        return true;
    else
        return false;
}

void Swap_Head::accept(){
    Move_Base::accept();
    path->put_worm_in_box(changed_particles,start_end);
    path->get_beads()->remove_old_list_references(rows_to_remove);
}

void Swap_Head::reject(){
     // std::cout << "Rejected" << std::endl;
    
    if(m != 0){
        
        std::vector<std::pair<int, int> > ht = path->get_beads()->get_old_worm_indices();
        int old_head_col_dist = ht[0].second-m+1;
        int old_tail_row = ht[1].first;
        if(old_head_col_dist < 0){
            old_head_col_dist=old_head_col_dist+path->get_parameters()->get_num_timeslices();
            old_tail_row++;
        }
        
        while((path->get_beads()->get_worm_indices()[1].first != old_tail_row)||(path->get_beads()->get_worm_indices()[0].second != old_head_col_dist))
            path->get_beads()->worm_pop_front();
        while((path->get_beads()->get_worm_indices()[1].first != ht[1].first)||(path->get_beads()->get_worm_indices()[0].second != ht[0].second))
            path->get_beads()->worm_pop_front(true);
    }
    
    Move_Base::reject();
    
}

/*******************************************************************
 
 SWAP TAIL MOVE
 
 ******************************************************************/

Swap_Tail::Swap_Tail(boost::shared_ptr<Path> path) : Move_Base(path){
    worm_move = true;
    worm_nec = true;
    
    move_name = "Swap Tail";
}

void Swap_Tail::attempt(){
    
    Move_Base::attempt();
    
    if(num_particles == 0){
        m = 0;
        reject();
        return;
    }
    
    start_end.resize(0);
    
    // std::cout << "Swap Tail Attempt" << std::endl;
    
    std::vector<std::pair<int, int> > ht = path->get_beads()->get_worm_indices();
    int tail_row = ht[1].first;
    int tail_col = ht[1].second;
    
    m = path->get_parameters()->get_mbar();
    
    
    std::vector<std::pair<size_t, dVector> > L_I = path->get_beads()->get_worm_tail_neighbors(path->get_parameters()->get_mbar());
    std::vector<std::pair<size_t, double> > rho0s_I;
    sig_I = 0;
    for(std::vector<std::pair<size_t, dVector> >::iterator it = L_I.begin(); it != L_I.end(); it++){
        ddVector pair;
        pair.push_back(path->get_beads()->get_worm_bead_data(tail_row, tail_col));
        pair.push_back(it->second);
        double rho0 = 1.0/pow((4*M_PI*path->get_parameters()->get_lambda()*path->get_parameters()->get_mbar()*path->get_parameters()->get_tau()),path->get_parameters()->get_ndim()/2.)*exp(-ka->get_action_pos(pair, path->get_parameters()->get_mbar()));
        rho0s_I.push_back(std::pair<size_t, double>(it->first, rho0));
        sig_I += rho0;
    }
    
    double rn = path->get_util()->randnormed(sig_I);
    size_t choice = 0;
    double accumulate = 0;
    for(std::vector<std::pair<size_t, double> >::iterator it = rho0s_I.begin(); it != rho0s_I.end(); it++){
        accumulate += it->second;
        if(rn < accumulate){
            choice = it->first;
            break;
        }
    }
    
    if(choice == 0){
        m = 0;
        reject();
        return;
    }
    
    ptcl = path->get_beads()->get_bead_indices(choice).first;
    changed_particles.push_back(ptcl);
    
    if(tail_col + path->get_parameters()->get_mbar() > path->get_parameters()->get_num_timeslices()){
        changed_particles = path->get_beads()->get_pair_rows(ptcl, tail_col, path->get_parameters()->get_mbar());
    }
    std::sort(changed_particles.begin(), changed_particles.end());
    auto last = std::unique(changed_particles.begin(), changed_particles.end());
    changed_particles.erase(last, changed_particles.end());
    
    for(int a = 0; a < path->get_parameters()->get_mbar(); a++){
        int slice = (tail_col + a)%path->get_parameters()->get_num_timeslices();
        if(changed_particles.size() == 1)
            old_action += pa->get_action_single_particle(ptcl, slice, true);
        else
            old_action += pa->get_action_multiple_particles(changed_particles, slice, true);
        old_action += pa->get_action(slice,0,true);
    }
    
    
    changed_particles.resize(0);
    
    size_t ksi_key = path->get_beads()->get_prev_bead_key(choice, m);
    if(!path->get_beads()->check_if_neighbor(path->get_beads()->get_worm_key(tail_row, tail_col), ksi_key)){
        m = 0;
        reject();
        return;
    }
    
    std::vector<std::pair<size_t, dVector> > L_ksi = path->get_beads()->get_neighbors(ksi_key, m);
    
    sig_ksi = 0;
    for(std::vector<std::pair<size_t, dVector> >::iterator it = L_ksi.begin(); it != L_ksi.end(); it++){
        ddVector pair;
        pair.push_back(path->get_beads()->get_bead_data(path->get_beads()->get_bead_indices(ksi_key).first, path->get_beads()->get_bead_indices(ksi_key).second));
        pair.push_back(it->second);
        double rho0 = 1.0/pow((4*M_PI*path->get_parameters()->get_lambda()*path->get_parameters()->get_mbar()*path->get_parameters()->get_tau()),path->get_parameters()->get_ndim()/2.)*exp(-ka->get_action_pos(pair, m));
        sig_ksi += rho0;
    }
    
    rows_to_remove = path->get_beads()->swap_into_worm_tail(path->get_beads()->get_bead_indices(choice).first, m);
    
    ddVector new_pos(m-1, dVector(path->get_parameters()->get_ndim(),0));
    dVector start_pos = path->get_beads()->get_worm_bead_data(ht[1].first, ht[1].second);
    int end_col = ht[1].second + m;
    int end_row = ht[1].first;
    if(end_col > path->get_parameters()->get_num_timeslices()){
        end_col = end_col%path->get_parameters()->get_num_timeslices();
        end_row++;
    }
    
    if(ht[1].first == end_row){
        changed_particles.push_back(end_row);
        start_end.push_back(std::pair<int,int>(ht[1].second+1, end_col-1));
    }
    else{
        if(ht[1].second != path->get_parameters()->get_num_timeslices()-1){
            changed_particles.push_back(ht[1].first);
            start_end.push_back(std::pair<int, int>(ht[1].second+1, path->get_parameters()->get_num_timeslices()-1));
        }
        if(end_col != 0){
            changed_particles.push_back(end_row);
            start_end.push_back(std::pair<int, int>(0, end_col-1));
        }
    }
    
    dVector end_pos = path->get_beads()->get_worm_bead_data(end_row, end_col);
    for(int n = 0; n < m-1; n++){
        double tau1 = (m-1-n)*path->get_parameters()->get_tau();
        dVector avex(path->get_parameters()->get_ndim());
        for(dVector::iterator it = avex.begin(); it!= avex.end(); it++){
            *it = (tau1*start_pos[it-avex.begin()]+path->get_parameters()->get_tau()*end_pos[it-avex.begin()])/(tau1+path->get_parameters()->get_tau());
            double box_size = path->get_parameters()->get_box_size();
            if(box_size != -1)
                *it = path->get_util()->per_bound_cond(*it, box_size);
        }
        double width = sqrt(2 * path->get_parameters()->get_lambda()/(1/path->get_parameters()->get_tau()+tau1));
        dVector bead_pos(path->get_parameters()->get_ndim());
        for(dVector::iterator it = bead_pos.begin(); it != bead_pos.end(); it++){
            *it = avex[it-bead_pos.begin()] + path->get_util()->randgaussian(width);
        }
        new_pos[n] = bead_pos;
        start_pos = bead_pos;
    }
    
    int counter = 0;
    int moving_tail_col = tail_col+1;
    int moving_tail_row = tail_row;
    while(counter < m-1){
        if(moving_tail_col >= path->get_beads()->get_worm_dims()[1]){
            moving_tail_row++;
            moving_tail_col = 0;
        }
        path->get_beads()->set_worm_bead_data(moving_tail_row, moving_tail_col, new_pos[counter]);
        moving_tail_col++;
        counter++;
    }
    
    for(int a = 0; a < path->get_parameters()->get_mbar(); a++){
        int slice = (tail_col + a)%path->get_parameters()->get_num_timeslices();
        new_action += pa->get_action(slice,0,true);
    }
    
    if(check_move()){        
        accept();
    }
    else
        reject();
    
}

bool Swap_Tail::check_move(){
    if(path->get_util()->randnormed(1)< sig_I/sig_ksi * exp(-(new_action-old_action)))
        return true;
    else
        return false;
}

void Swap_Tail::accept(){
    Move_Base::accept();
    path->put_worm_in_box(changed_particles, start_end);
    path->get_beads()->remove_old_list_references(rows_to_remove);

}

void Swap_Tail::reject(){
     // std::cout << "Rejected" << std::endl;
    if(m != 0){
        
        std::vector<std::pair<int, int> > ht = path->get_beads()->get_old_worm_indices();
        int old_tail_col_dist = ht[1].second+m-1;
        int old_tail_row = ht[1].first;
        if(old_tail_col_dist >=path->get_parameters()->get_num_timeslices()){
            old_tail_col_dist=old_tail_col_dist%path->get_parameters()->get_num_timeslices();
            old_tail_row++;
        }
        
        while((path->get_beads()->get_worm_indices()[1].first != old_tail_row)||(path->get_beads()->get_worm_indices()[1].second != old_tail_col_dist))
            path->get_beads()->worm_pop_back();
        while((path->get_beads()->get_worm_indices()[1].first != ht[1].first)||(path->get_beads()->get_worm_indices()[1].second != ht[1].second))
            path->get_beads()->worm_pop_back(true);
        
        Move_Base::reject();
    }
}
