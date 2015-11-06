//
//  Moves.cpp
//  PIMC
//
//  Created by Adith Ramamurti on 5/14/15.
//  Copyright (c) 2015 Adith Ramamurti. All rights reserved.
//

#include "moves.h"



bool moves::comMove(paths* path, int ptcl){
    float delta;
    if(path->getParam()->getBoxSize() == -1)
        delta = 0.5;
    else{
        delta = sqrt(path->getParam()->getlam()*path->getParam()->gettau());
    }
    
    vectorf shift(path->getParam()->getndim());
    for(vectorf::iterator it = shift.begin(); it != shift.end(); it++ )
        *it = path->getUte()->randgaussian(delta);
    
    path->getBeads()->setOld();
    
    float oldAct = 0.0;
    for(int slice = 0; slice < path->getParam()->getNumTimeSlices(); slice++){
        oldAct += path->potentialAction(slice);
    }
    
    
    path->getBeads()->shiftAll(ptcl, shift);
    
    float newAct = 0.0;
    for(int slice = 0; slice < path->getParam()->getNumTimeSlices(); slice++){
        newAct += path->potentialAction(slice);
    }
    
    if(path->getUte()->randnormed(1)<exp(-(newAct - oldAct))){
        path->putInBox();
        std::vector<int> lc = path->getLastChgd();
        lc.push_back(ptcl);
        path->setlastChgd(lc);
        return true;
    }
    else{
        path->getBeads()->resetOld();
        return false;
    }
}


inline void moves::bisectionMove(paths* path, int ptcl, int start, int m){
    if(m != 1 && m%2 == 0){
        int slice = (start + m/2);
        float tau1 = (m/2)*path->getParam()->gettau();
        vectorf move(0);
        vectorff bds = path->getBeads()->getPair(ptcl, start, m);
        vectorf aved = path->getUte()->avedist(bds,path->getParam()->getBoxSize());
        float width = sqrt(path->getParam()->getlam()*tau1);
        for(int ndim = 0; ndim < path->getParam()->getndim(); ndim++){
            move.push_back(aved[ndim] + path->getUte()->randgaussian(width));
        }
        path->getBeads()->setOne(ptcl, slice, move);
        bisectionMove(path, ptcl, start, m/2);
        bisectionMove(path, ptcl, slice, m/2);
    }
}

bool moves::bisectionMoveHelper(paths* path, int ptcl){
    int m = path->getDist();
    
    int start = path->getUte()->randint(path->getParam()->getNumTimeSlices());
    
    int ls = path->getlastLocs()[0];
    int le = path->getlastLocs()[1];
    
    std::vector<int> lc;
    
    if(le > path->getParam()->getNumTimeSlices() && start < le% path->getParam()->getNumTimeSlices())
        lc = path->getBeads()->getRSO();
    else
        lc = path->getLastChgd();
    
    if(path->getParam()->isboson())
        if((start+m > ls &&start+m <=le) || ((start >=ls) && (start <=le))|| (le > path->getParam()->getNumTimeSlices() && start < le % path->getParam()->getNumTimeSlices()))
            path->recompSingProb(lc, start);
    
    float oldPotAct = 0.0;
    
    
    path->getBeads()->setOld();
    
    
    for(int a = 0; a < m+1; a++){
        int slice = (start + a)%path->getParam()->getNumTimeSlices();
        oldPotAct += path->potentialAction(slice);
    }
    
    std::vector<int> identity(path->getParam()->getNumParticles());
    iota(identity.begin(),identity.end(),0);
    std::vector<int> chosenPerm = identity;
    std::vector<int> origpart(0);
    std::vector<int> permpart(0);
    
    if(path->getParam()->isboson()){
        chosenPerm = pickPermutation(path, start);
        
        for(std::vector<int>::iterator it = identity.begin(); it != identity.end(); it++)
            if(*it != chosenPerm[*it]){
                origpart.push_back(*it);
                permpart.push_back(chosenPerm[*it]);
            }
    }
    
    path->getBeads()->setswap(origpart, permpart, start, m);
    path->getBeads()->swap();
    
    if(permpart.size() == 0){
        bisectionMove(path, ptcl, start, m);
        permpart.push_back(ptcl);
    }
    else
        for(std::vector<int>::iterator ptcl = origpart.begin(); ptcl !=origpart.end(); ptcl++)
            bisectionMove(path, *ptcl, start, m);
    
    float newPotAct = 0.0;
    for(int a = 0; a < m+1; a++ ){
        int slice = (start + a)%path->getParam()->getNumTimeSlices();
        newPotAct += path->potentialAction(slice);
    }
    
    float potDiff = newPotAct-oldPotAct;
    float rn = path->getUte()->randnormed(1);
    
    
    if(exp(-potDiff) < rn){
        
        path->getBeads()->swap(true);
        path->getBeads()->resetOld();
        path->getBeads()->resetSwap();
        
        return false;
        
    }
        
    path->setlastChgd(permpart);
    path->slstep(start, start+m);
    path->getBeads()->setPrevSwap();
    path->putInBox();
    return true;
}

inline std::vector<int> moves::pickPermutation(paths* path, int start){
    vectorf permWeight = (*path->getProbList())[start];
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
    
    std::vector<int> chosenPerm = (*path->getPermList())[choice];
    
    return chosenPerm;
}
