//
//  Moves.cpp
//  PIMC
//
//  Created by Adith Ramamurti on 5/14/15.
//  Copyright (c) 2015 Adith Ramamurti. All rights reserved.
//

#include "moves.h"



bool moves::comMove(Path* path, int ptcl){
    float delta;
    if(path->get_parameters()->get_box_size() == -1)
        delta = 0.5;
    else{
        delta = sqrt(path->get_parameters()->get_lambda()*path->get_parameters()->get_tau());
    }
    
    vectorf shift(path->get_parameters()->get_ndim());
    for(vectorf::iterator it = shift.begin(); it != shift.end(); it++ )
        *it = path->getUte()->randgaussian(delta);
    
    path->get_beads()->setOld();
    
    float oldAct = 0.0;
    for(int slice = 0; slice < path->get_parameters()->get_num_timeslices(); slice++){
        oldAct += path->potentialAction(slice);
    }
    
    
    path->get_beads()->shiftAll(ptcl, shift);
    
    float newAct = 0.0;
    for(int slice = 0; slice < path->get_parameters()->get_num_timeslices(); slice++){
        newAct += path->potentialAction(slice);
    }
    
    if(path->getUte()->randnormed(1)<exp(-(newAct - oldAct))){
        path->put_in_box();
        std::vector<int> lc = path->get_last_changed();
        lc.push_back(ptcl);
        path->set_last_changed(lc);
        return true;
    }
    else{
        path->get_beads()->resetOld();
        return false;
    }
}


inline void moves::bisectionMove(Path* path, int ptcl, int start, int m){
    if(m != 1 && m%2 == 0){
        int slice = (start + m/2);
        float tau1 = (m/2)*path->get_parameters()->get_tau();
        vectorf move(0);
        vectorff bds = path->get_beads()->getPair(ptcl, start, m);
        vectorf aved = path->getUte()->avedist(bds,path->get_parameters()->get_box_size());
        float width = sqrt(path->get_parameters()->get_lambda()*tau1);
        for(int ndim = 0; ndim < path->get_parameters()->get_ndim(); ndim++){
            move.push_back(aved[ndim] + path->getUte()->randgaussian(width));
        }
        path->get_beads()->setOne(ptcl, slice, move);
        bisectionMove(path, ptcl, start, m/2);
        bisectionMove(path, ptcl, slice, m/2);
    }
}

bool moves::bisectionMoveHelper(Path* path, int ptcl){
    int m = path->getDist();
    
    int start = path->getUte()->randint(path->get_parameters()->get_num_timeslices());
    
    int ls = path->get_last_locs()[0];
    int le = path->get_last_locs()[1];
    
    std::vector<int> lc;
    
    if(le > path->get_parameters()->get_num_timeslices() && start < le% path->get_parameters()->get_num_timeslices())
        lc = path->get_beads()->getRSO();
    else
        lc = path->get_last_changed();
    
    if(path->get_parameters()->is_boson())
        if((start+m > ls &&start+m <=le) || ((start >=ls) && (start <=le))|| (le > path->get_parameters()->get_num_timeslices() && start < le % path->get_parameters()->get_num_timeslices()))
            path->slice_perm_prob(lc, start);
    
    float oldPotAct = 0.0;
    
    
    path->get_beads()->setOld();
    
    
    for(int a = 0; a < m+1; a++){
        int slice = (start + a)%path->get_parameters()->get_num_timeslices();
        oldPotAct += path->potentialAction(slice);
    }
    
    std::vector<int> identity(path->get_parameters()->get_num_particles());
    iota(identity.begin(),identity.end(),0);
    std::vector<int> chosenPerm = identity;
    std::vector<int> origpart(0);
    std::vector<int> permpart(0);
    
    if(path->get_parameters()->is_boson()){
        chosenPerm = pickPermutation(path, start);
        
        for(std::vector<int>::iterator it = identity.begin(); it != identity.end(); it++)
            if(*it != chosenPerm[*it]){
                origpart.push_back(*it);
                permpart.push_back(chosenPerm[*it]);
            }
    }
    
    path->get_beads()->setswap(origpart, permpart, start, m);
    path->get_beads()->swap();
    
    if(permpart.size() == 0){
        bisectionMove(path, ptcl, start, m);
        permpart.push_back(ptcl);
    }
    else
        for(std::vector<int>::iterator ptcl = origpart.begin(); ptcl !=origpart.end(); ptcl++)
            bisectionMove(path, *ptcl, start, m);
    
    float newPotAct = 0.0;
    for(int a = 0; a < m+1; a++ ){
        int slice = (start + a)%path->get_parameters()->get_num_timeslices();
        newPotAct += path->potentialAction(slice);
    }
    
    float potDiff = newPotAct-oldPotAct;
    float rn = path->getUte()->randnormed(1);
    
    
    if(exp(-potDiff) < rn){
        
        path->get_beads()->swap(true);
        path->get_beads()->resetOld();
        path->get_beads()->resetSwap();
        
        return false;
        
    }
        
    path->set_last_changed(permpart);
    path->set_last_step(start, start+m);
    path->get_beads()->setPrevSwap();
    path->put_in_box();
    return true;
}

inline vectori moves::pickPermutation(Path* path, int start){
    vectorf permWeight = (*path->get_prob_list())[start];
    vectorf::iterator it2;
    
    float sum = 0.0;
    for(it2 = permWeight.begin(); it2 != permWeight.end(); it2++)
        sum += *it2;
    
    int choice = 0;
    float rn = path->getUte()->randnormed(sum);
    float probsum = 0.0;
    for(int i = 0; i < permWeight.size(); i++){
        probsum += permWeight[i];
        if(rn<=probsum){
            choice = i;
            break;
        }
    }
    
    std::vector<int> chosenPerm = (*path->get_perm_list())[choice];
    
    return chosenPerm;
}
