//
//  Moves.cpp
//  PIMC
//
//  Created by Adith Ramamurti on 5/14/15.
//  Copyright (c) 2015 Adith Ramamurti. All rights reserved.
//

#include "moves.h"



bool moves::comMove(paths* path, int ptcl){
    double delta;
    if(path->getParam()->getBoxSize() == -1)
        delta = 0.5;
    else
        delta = path->getParam()->getBoxSize()/10.0;
    
    std::vector<double> shift(path->getParam()->getndim());
    for(std::vector<double>::iterator it = shift.begin(); it != shift.end(); it++ )
        *it = delta *(path->getUte()->randnormed(2)-1);
    
    paths::list_ptr oldLL(new LinkedList<std::vector<double>>(*path->getBeads()));
    
    double oldAct = 0.0;
    for(int slice = 0; slice < path->getParam()->getNumTimeSlices(); slice++){
        oldAct += path->potentialAction(slice);
    }
    
    
    path->getBeads()->shiftAll(ptcl, shift);
    
    double newAct = 0.0;
    for(int slice = 0; slice < path->getParam()->getNumTimeSlices(); slice++){
        newAct += path->potentialAction(slice);
    }
    
    if(path->getUte()->randnormed(1)<exp(-(newAct - oldAct))){
        return true;
    }
    else{
        path->setBeads(*oldLL);
        return false;
    }
}

bool moves::stagingMoveHelper(paths* path, int ptcl){
    int m = path->getDist();
    
    int start = path->getUte()->randint(path->getParam()->getNumTimeSlices());
    
    
    double permTot0 = path->recompSingProb(start);
    
    double oldPotAct = 0.0;
    
    paths::list_ptr oldLL(new LinkedList<std::vector<double>>(*path->getBeads()));
    
    for(int a = 0; a < m+1; a++){
        int slice = (start + a)%path->getParam()->getNumTimeSlices();
        oldPotAct += path->potentialAction(slice);
    }
    
    std::vector<int> identity(path->getParam()->getNumParticles());
    iota(identity.begin(),identity.end(),0);
    std::vector<int> chosenPerm = identity;
    std::vector<int> origpart(0);
    std::vector<int> permpart(0);
    
    if(path->getParam()->getboson()){
        chosenPerm = pickPermutation(path, start);
        //std::cout << start << "\t" << start+m<< std::endl;;
        
        for(std::vector<int>::iterator it = identity.begin(); it != identity.end(); it++)
            if(*it != chosenPerm[*it]){
                origpart.push_back(*it);
                permpart.push_back(chosenPerm[*it]);
                //std::cout << *it << "\t" << chosenPerm[*it]<< std::endl;;
            }
    }
    
    path->getBeads()->swap(origpart, permpart, start, m);
    
    if(permpart.size() == 0)
        stagingMove(path, ptcl, start, m);
    else
        for(std::vector<int>::iterator ptcl = origpart.begin(); ptcl !=origpart.end(); ptcl++)
            stagingMove(path, *ptcl, start, m);
    
    double newPotAct = 0.0;
    for(int a = 0; a < m+1; a++ ){
        int slice = (start + a)%path->getParam()->getNumTimeSlices();
        newPotAct += path->potentialAction(slice);
    }
    
    double potDiff = newPotAct-oldPotAct;
    
    if(exp(-potDiff) < path->getUte()->randnormed(1)){
        //std::cout<<"rejected"<<std::endl;
        
        path->getBeads()->reverseswap(origpart, permpart, start, m);
        path->setBeads(*oldLL);
        return false;
    }
    
    double permTot1 = path->recompSingProb(start);
    
    if((permTot0/permTot1) < path->getUte()->randnormed(1)){
        path->getBeads()->reverseswap(origpart, permpart, start, m);
        path->setBeads(*oldLL);
        return false;
    }
    return true;
}


void moves::stagingMove(paths *path, int ptcl, int start, int m){
    int end = (start + m)%path->getParam()->getNumTimeSlices();

    for(int a = 1; a < m; a++){
        int slice = (start + a)%path->getParam()->getNumTimeSlices();
        int slicem1 = (slice - 1 + path->getParam()->getNumTimeSlices())%path->getParam()->getNumTimeSlices();
        double tau1 = (m-a)*path->getParam()->gettau();
        
        std::vector<double> move(0);
        for(int ndim = 0; ndim < path->getParam()->getndim();ndim++){
            double avex = (tau1*path->getBeads()->getOne(ptcl, slicem1)[ndim]+path->getParam()->gettau()*path->getBeads()->getOne(ptcl, end)[ndim])/(path->getParam()->gettau()+tau1);
            double width = sqrt(2.0*path->getParam()->getlam()/(1.0/path->getParam()->gettau() + 1.0/tau1));
            move.push_back(avex + path->getUte()->randgaussian(width));
        }
        path->getBeads()->setOne(ptcl, slice, move);
    }
    
}

void moves::bisectionMove(paths* path, int ptcl, int start, int m){
    int end = (start + m)%path->getParam()->getNumTimeSlices();

    if(m != 1 && m%2 == 0){
        int slice = (start + m/2)%path->getParam()->getNumTimeSlices();
        double tau1 = (m/2)*path->getParam()->gettau();
        std::vector<double> move(0);
        for(int ndim = 0; ndim < path->getParam()->getndim(); ndim++){
            double avex = (path->getBeads()->getOne(ptcl, start)[ndim]+path->getBeads()->getOne(ptcl, end)[ndim])/(2);
            double width = sqrt(path->getParam()->getlam()*tau1);
            move.push_back(avex + path->getUte()->randgaussian(width));
        }
        path->getBeads()->setOne(ptcl, slice, move);
        bisectionMove(path, ptcl, start, m/2);
        bisectionMove(path, ptcl, slice, m/2);
    }
}

bool moves::bisectionMoveHelper(paths* path, int ptcl){
    int m = path->getDist();
    
    int start = path->getUte()->randint(path->getParam()->getNumTimeSlices());
    
    double permTot0 = path->recompSingProb(start);
    
    double oldPotAct = 0.0;
    
    paths::list_ptr oldLL(new LinkedList<std::vector<double>>(*path->getBeads()));
    
    for(int a = 0; a < m+1; a++){
        int slice = (start + a)%path->getParam()->getNumTimeSlices();
        oldPotAct += path->potentialAction(slice);
    }
    
    std::vector<int> identity(path->getParam()->getNumParticles());
    iota(identity.begin(),identity.end(),0);
    std::vector<int> chosenPerm = identity;
    std::vector<int> origpart(0);
    std::vector<int> permpart(0);
    
    if(path->getParam()->getboson()){
        chosenPerm = pickPermutation(path, start);
        //std::cout << start << "\t" << start+m<< std::endl;;
        
        for(std::vector<int>::iterator it = identity.begin(); it != identity.end(); it++)
            if(*it != chosenPerm[*it]){
                origpart.push_back(*it);
                permpart.push_back(chosenPerm[*it]);
                //std::cout << *it << "\t" << chosenPerm[*it]<< std::endl;;
            }
    }
    
    path->getBeads()->swap(origpart, permpart, start, m);
    
    
    if(permpart.size() == 0)
        bisectionMove(path, ptcl, start, m);
    else
        for(std::vector<int>::iterator ptcl = origpart.begin(); ptcl !=origpart.end(); ptcl++)
            bisectionMove(path, *ptcl, start, m);
    
    double newPotAct = 0.0;
    for(int a = 0; a < m+1; a++ ){
        int slice = (start + a)%path->getParam()->getNumTimeSlices();
        newPotAct += path->potentialAction(slice);
    }
    
    double potDiff = newPotAct-oldPotAct;
    
    if(exp(-potDiff) < path->getUte()->randnormed(1)){
        //std::cout<<"rejected"<<std::endl;
        
        path->getBeads()->reverseswap(origpart, permpart, start, m);
        path->setBeads(*oldLL);
        return false;
    }
    
    double permTot1 = path->recompSingProb(start);
    
    if((permTot0/permTot1) < path->getUte()->randnormed(1)){
        path->getBeads()->reverseswap(origpart, permpart, start, m);
        path->setBeads(*oldLL);
        return false;
    }
    return true;
}

std::vector<int> moves::pickPermutation(paths* path, int start){
    std::vector<double> permWeight = (*path->getProbList())[start];
    std::vector<double>::iterator it2;
    
    double sum = 0.0;
    for(it2 = permWeight.begin(); it2 != permWeight.end(); it2++)
        sum += *it2;
    
    int choice = 0;
    double rn = path->getUte()->randnormed(sum);
    double probsum = 0.0;
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
