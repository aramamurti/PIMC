//
//  pimc.cpp
//  PIMC
//
//  Created by Adith Ramamurti on 5/20/15.
//  Copyright (c) 2015 Adith Ramamurti. All rights reserved.
//

#include "pimc.h"
#include "moves.h"
#include <iostream>

using namespace std;

pimc::pimc(){
    skip = 100;
    equil = 1000;
    numacceptc = 0;
    numaccepts = 0;
}

vector<double> pimc::run(int numSteps, paths* path, vector<bool> pmv){
    vector<double> energytr = {};
    
    moves mvs;
    
    for(int step = 0; step < numSteps; step++){
        int ptcl = (int) path->getUte()->randnormed(path->getParam()->getNumParticles())%path->getParam()->getNumParticles();
        if(pmv[0])
            if(mvs.comMove(path, ptcl))
                numacceptc += 1;
        if(pmv[1])
            if(mvs.stagingMove(path, ptcl))
                numaccepts += 1;
        if(pmv[2])
            if(mvs.bisectionMoveHelper(path, ptcl))
                numaccepts += 1;
        
        if(step % skip == 0 && step>equil)
            energytr.push_back(path->energy());
    }
    
    cout << "\nCenter of mass acceptance: "<< 1.0*numacceptc/(numSteps*path->getParam()->getNumParticles()) << "\n" << "Staging acceptance: "<< 1.0*numaccepts/(numSteps*path->getParam()->getNumParticles()) << "\n\n";
    
    return energytr;
}