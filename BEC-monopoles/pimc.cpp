//
//  pimc.cpp
//  PIMC
//
//  Created by Adith Ramamurti on 5/20/15.
//  Copyright (c) 2015 Adith Ramamurti. All rights reserved.
//

#include "moves.h"
#include "paths.h"
#include "pimc.h"
#include "uni_header.h"
#include <unistd.h>



pimc::pimc(){
    numacceptc = 0;
    numaccepts = 0;
    numacceptb = 0;
}

std::vector<double> pimc::run(int numSteps, paths* path, std::ofstream &f){
    std::vector<double> energytr(0);
    
    moves mvs;
    
    for(int step = 0; step < numSteps; step++){
        int ptcl = (int) path->getUte()->randnormed(path->getParam()->getNumParticles())%path->getParam()->getNumParticles();
        if(mvs.comMove(path, ptcl))
            numacceptc += 1;
        /*if(mvs.stagingMoveHelper(path, ptcl))
            numaccepts += 1;*/
        if(mvs.stagingMoveHelper(path, ptcl))
            numacceptb += 1;
        
        if(step % path->getParam()->getSkip() == 0 && step>path->getParam()->getEquil()){
            double en = path->energy();
            f << step << "\t" << en << "\n";
            energytr.push_back(en);
        }
    }
    
    f << "\nCenter of mass acceptance: "<< 1.0*numacceptc/(numSteps) << "\nStaging acceptance: "<< 1.0*numaccepts/(numSteps)<< "\nBisection acceptance: "<< 1.0*numacceptb/(numSteps) << "\n\n";
    
    return energytr;
}