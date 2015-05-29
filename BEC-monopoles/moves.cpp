//
//  Moves.cpp
//  PIMC
//
//  Created by Adith Ramamurti on 5/14/15.
//  Copyright (c) 2015 Adith Ramamurti. All rights reserved.
//

#include "moves.h"
#include <vector>
#include <cmath>
#include <iostream>
#include <algorithm>

using namespace std;

bool moves::comMove(paths* path, int ptcl){
    double delta = 0.5;
    double shift[] = {delta *(path->getUte()->randnormed(2)-1),delta *(path->getUte()->randnormed(2)-1),delta *(path->getUte()->randnormed(2)-1)};
    
    vector<vector<double>> bptc = {};
    for(int i = 0; i<path->beads.size();i++){
        bptc.push_back(path->beads[i][ptcl]);
    }
    
    double oldAct = 0.0;
    for(int slice = 0; slice < path->getParam()->getNumTimeSlices(); slice++){
        oldAct += path->potentialAction(slice);
    }
    for(int slice = 0; slice < path->getParam()->getNumTimeSlices(); slice++){
        for(int ndim = 0; ndim < path->getParam()->getndim(); ndim++)
            path->beads[slice][ptcl][ndim] = bptc[slice][ndim]+shift[ndim];
    }
    double newAct = 0.0;
    for(int slice = 0; slice < path->getParam()->getNumTimeSlices(); slice++){
        newAct += path->potentialAction(slice);
    }
    
    if(path->getUte()->randnormed(1)<exp(-(newAct - oldAct))){
        return true;
    }
    else{
        for(int i = 0; i<path->beads.size();i++){
            path->beads[i][ptcl] = bptc[i];
        }
        return false;
    }
}

bool moves::stagingMove(paths* path, int ptcl){
    int m = 16;
    
    int alpha_start = path->getUte()->randint(path->getParam()->getNumTimeSlices());
    int alpha_end = (alpha_start + m)%path->getParam()->getNumTimeSlices();
    
    vector<vector<double>> oldbeads(m-1, vector<double>(path->getParam()->getndim(),0.0));
    double oldAct = 0.0;
    for(int a = 1; a <m; a++){
        int slice = (alpha_start + a)%path->getParam()->getNumTimeSlices();
        oldbeads[a-1]= path->beads[slice][ptcl];
        oldAct += path->potentialAction(slice);
    }
    
    if(path->getParam()->getboson())
        pickPermutation(path, alpha_start, alpha_end);
    
    double newAct = 0.0;
    for(int a = 1; a <m; a++){
        int slice = (alpha_start + a)%path->getParam()->getNumTimeSlices();
        int slicem1 = (slice + path->getParam()->getNumTimeSlices() - 1)%path->getParam()->getNumTimeSlices();
        double tau1 = (m-a)*path->getParam()->gettau();
        for(int ndim = 0; ndim < path->getParam()->getndim();ndim++){
            double avex = (tau1*path->beads[slicem1][ptcl][ndim]+path->getParam()->gettau()*path->beads[alpha_end][ptcl][ndim])/(path->getParam()->gettau()+tau1);
            double width = sqrt(2.0*path->getParam()->getlam()/(1.0/path->getParam()->gettau() + 1.0/tau1));
            double move = avex + path->getUte()->randgaussian(width);
            path->beads[slice][ptcl][ndim] = move;
        }
        newAct += path->potentialAction(slice);
    }
    
    if(path->getUte()->randnormed(1)<exp(-(newAct - oldAct))){
        return true;
    }
    else{
        for(int a = 1; a< m; a++){
            int slice = (alpha_start + a)%path->getParam()->getNumTimeSlices();
            path->beads[slice][ptcl] = oldbeads[a-1];
        }
        return false;
    }
}

void moves::bisectionMove(paths* path, int ptcl, int alpha_start, int alpha_end, int m){
    
    if(m != 1 && m%2 == 0){
        int slice = (alpha_start + m/2)%path->getParam()->getNumTimeSlices();
        double tau1 = (m/2)*path->getParam()->gettau();
        for(int ndim = 0; ndim < path->getParam()->getndim(); ndim++){
            double avex = (path->beads[alpha_start][ptcl][ndim]+path->beads[alpha_end][ptcl][ndim])/2;
            double width = sqrt(path->getParam()->gettau()*tau1);
            double move = avex + path->getUte()->randgaussian(width);
            path->beads[slice][ptcl][ndim] = move;
        }
        bisectionMove(path, ptcl, alpha_start, slice, m/2);
        bisectionMove(path, ptcl, slice, alpha_end, m/2);
    }
}

bool moves::bisectionMoveHelper(paths* path, int ptcl){
    int m = 16;

    int alpha_start = path->getUte()->randint(path->getParam()->getNumTimeSlices());
    int alpha_end = (alpha_start + m)%path->getParam()->getNumTimeSlices();
    
    vector<vector<double>> oldbeads(m-1, vector<double>(path->getParam()->getndim(),0.0));
    
    double oldKinAct = 0.0;
    double oldPotAct = 0.0;
    for(int a = 1; a < m; a++){
        int slice = (alpha_start + a)%path->getParam()->getNumTimeSlices();
        oldbeads[a-1]= path->beads[slice][ptcl];
        oldPotAct += path->potentialAction(slice);
        oldKinAct += path->kineticAction(slice,1);
    }
    vector<int> permpart;
    
    if(path->getParam()->getboson())
        permpart = pickPermutation(path, alpha_start, alpha_end);
    
    double newPotAct = 0.0;

    if(permpart.size() == 0)
        bisectionMove(path, ptcl, alpha_start, alpha_end, m);
    else
        for(int part = 0; part < permpart.size(); part++)
            bisectionMove(path, permpart[part], alpha_start, alpha_end, m);
    
    double newKinAct = 0.0;
    for(int a = 1; a <= m; a++ ){
        int slice = (alpha_start + a)%path->getParam()->getNumTimeSlices();
        newPotAct += path->potentialAction(slice);
        newKinAct += path->kineticAction(slice, 1);
    }
    
    double potDiff = newPotAct-oldPotAct;
    double kinDiff = newKinAct-oldKinAct;
    
    double diff = 0;
    if(permpart.size() == 0)
        diff = potDiff;
    else{
        diff = potDiff+kinDiff;
    }
    
    if(path->getUte()->randnormed(1)<exp(-diff))
        return true;
    else{
        for(int a = 1; a< m; a++){
            int slice = (alpha_start + a)%path->getParam()->getNumTimeSlices();
            path->beads[slice][ptcl] = oldbeads[a-1];
        }
        return false;
    }
}

vector<int> moves::pickPermutation(paths* path, int alpha_start, int alpha_end){
    int dist = (alpha_end-alpha_start+path->getParam()->getNumTimeSlices())%path->getParam()->getNumTimeSlices();
    vector<double> permWeight;
    vector<double> problist;
    
    vector<vector<int>> permVec(path->getPermList());
    
    vector<vector<double>> oldEndBeads(path->beads[alpha_end]);

    for(int i = 0; i < permVec.size(); i++){
        vector<int> oneperm = permVec[i];
        for(int ptcl = 0; ptcl < path->getParam()->getNumParticles(); ptcl++){
            path->beads[alpha_end][ptcl] = oldEndBeads[oneperm[ptcl]];
        }
        
        double action = 0.0;
        for(int ptcl = 0; ptcl < path->getParam()->getNumParticles(); ptcl++){
            action += path->kineticAction(alpha_start,dist);
        }
        problist.push_back(action);
        path->beads[alpha_end] = oldEndBeads;
    }
    
    double sum = 0.0;
    vector<double>::iterator it2;
    for(it2 = problist.begin(); it2 != problist.end(); it2++)
        sum += *it2;
    
    for(it2 = problist.begin(); it2 != problist.end(); it2++){
        *it2 /= sum;
        permWeight.push_back(*it2);
    }
    
    double checkSum = 0.0;
    for(int i = 0; i < permWeight.size(); i++)
        checkSum += permWeight[i];
    
    int choice = 0;
    if(checkSum == 1){
        double rn = path->getUte()->randnormed(1);
        double probsum = 0.0;
        for(int i = 0; i < permWeight.size(); i++){
            probsum += permWeight[i];
            if(rn<=probsum)
                choice = i;
        }
    }
        
    vector<int> chosenPerm = permVec[choice];
    
    vector<int> chgdpart = {};
    for(int slice = alpha_end; slice < path->getParam()->getNumTimeSlices(); slice ++){
        vector<vector<double>> beforePerm(path->beads[slice]);
        for(int ptclnum = 0; ptclnum < path->getParam()->getNumParticles(); ptclnum++){
            path->beads[slice][ptclnum] = beforePerm[chosenPerm[ptclnum]];
            if(ptclnum != chosenPerm[ptclnum]){
                chgdpart.push_back(ptclnum);
            }
        }
    }
    return chgdpart;
}
