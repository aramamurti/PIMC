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
    
    /*std::ofstream p1;
    std::ofstream p2;
    std::ofstream p3;
    std::ofstream p4;
    
    p1.open("p1.csv");
    p2.open("p2.csv");
    p3.open("p3.csv");
    p4.open("p4.csv");*/
    
    f << "Step No." << "\t\t" << "Energy/atom" << "\t\t" << "KE/atom" << "\t\t" << "PE/atom" << std::endl;

    
    for(int step = 0; step < numSteps; step++){
        int ptcl = (int) path->getUte()->randnormed(path->getParam()->getNumParticles())%path->getParam()->getNumParticles();
        if(mvs.comMove(path, ptcl))
            numacceptc += 1;
        /*if(mvs.stagingMoveHelper(path, ptcl))
            numaccepts += 1;*/
        if(mvs.bisectionMoveHelper(path, ptcl))
            numacceptb += 1;
        if(step % path->getParam()->getSkip() == 0 && step >= path->getParam()->getEquil()){
            double en = path->energy();

            /*if(step == path->getParam()->getEquil()){
                path->getBeads()->printLinkedList(0, p1);
                path->getBeads()->printLinkedList(1, p2);
                path->getBeads()->printLinkedList(2, p3);
                path->getBeads()->printLinkedList(3, p4);
                p1.close();
                p2.close();
                p3.close();
                p4.close();
            }*/

            f << step << "\t\t" << en/path->getParam()->getNumParticles() << "\t\t" << path->kineticEnergy()/path->getParam()->getNumParticles() << "\t\t" << path->potentialEnergy()/path->getParam()->getNumParticles() << std::endl;
            energytr.push_back(en);
        }
    }
    
    f << "\nCenter of mass acceptance: "<< 1.0*numacceptc/(numSteps) << "\nStaging acceptance: "<< 1.0*numaccepts/(numSteps)<< "\nBisection acceptance: "<< 1.0*numacceptb/(numSteps) << "\n" << std::endl;
    
    return energytr;
}