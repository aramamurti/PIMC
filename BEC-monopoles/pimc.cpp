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

vector<double> pimc::run(int numSteps, paths* path){
    vector<double> energytr = {};
    
    moves mvs;
        
    for(int step = 0; step < numSteps; step++){
        int ptcl = (int) path->getUte()->randnormed(path->getNumParticles())%path->getNumParticles();
        if(mvs.comMove(path, ptcl))
            numacceptc += 1;
        
        /*if(mvs.stagingMove(path, ptcl))
            numaccepts += 1;*/
        
        if(mvs.bisectionMoveHelper(path, ptcl))
            numaccepts += 1;
        
        if(step % skip == 0 && step>equil)
            energytr.push_back(path->energy());
    }
    
    cout << "Center of mass acceptance: "<< 1.0*numacceptc/(numSteps*path->getNumParticles()) << "\n";
    cout << "Staging acceptance: "<< 1.0*numaccepts/(numSteps*path->getNumParticles()) << "\n\n";
    
    return energytr;
}