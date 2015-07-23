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

void pimc::run(int numSteps, paths* path, std::ofstream &f1, std::ofstream &f2, std::vector<double> &energytr, std::vector<std::vector<int>> &cycleList){
    
    moves mvs;
        
    f1 << "Step No." << ", " << "Energy/atom" << ", " << "KE/atom" << ", " << "PE/atom" << std::endl;

    
    for(int step = 0; step < numSteps; step++){
        int ptcl = (int) path->getUte()->randnormed(path->getParam()->getNumParticles())%path->getParam()->getNumParticles();
        if(mvs.comMove(path, ptcl))
            numacceptc += 1;
        if(mvs.bisectionMoveHelper(path, ptcl))
            numacceptb += 1;
        if(step % path->getParam()->getSkip() == 0 && step >= path->getParam()->getEquil()){
            double en = path->energy();
            f1 << step << ", " << en/path->getParam()->getNumParticles() << ", " << path->kineticEnergy()/path->getParam()->getNumParticles() << ", " << path->potentialEnergy()/path->getParam()->getNumParticles() << std::endl;
            std::vector<int> cycles = path->getCycles();
            for(std::vector<int>::iterator it = cycles.begin(); it != cycles.end(); it++){
                f2 << *it;
                if(cycles.size() - (it-cycles.begin()) != 1)
                    f2 << ", ";
            }
            f2 << std::endl;
            energytr.push_back(en);
            cycleList.push_back(cycles);
        }
    }
    
    f1 << "\nCenter of mass acceptance: "<< 1.0*numacceptc/(numSteps) << "\nBisection acceptance: "<< 1.0*numacceptb/(numSteps) << "\n" << std::endl;
    
}