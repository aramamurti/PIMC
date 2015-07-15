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

std::vector<double> pimc::run(int numSteps, paths* path, Gnuplot &g){
    std::vector<double> energytr(0);
    
    moves mvs;
    /*g.set_style("lines");
    
    Gnuplot g1("plot2");
    g1.set_style("lines");*/
    
    for(int step = 0; step < numSteps; step++){
        std::cout<<step<<std::endl;
        int ptcl = (int) path->getUte()->randnormed(path->getParam()->getNumParticles())%path->getParam()->getNumParticles();
        if(mvs.comMove(path, ptcl))
            numacceptc += 1;
        /*if(mvs.stagingMoveHelper(path, ptcl))
            numaccepts += 1;*/
        if(mvs.stagingMoveHelper(path, ptcl))
            numacceptb += 1;
        
        if(step % path->getParam()->getSkip() == 0 && step>path->getParam()->getEquil()){
            energytr.push_back(path->energy());
            /*g.reset_plot();
            g.plot_x(energytr, "energy");
            
            g1.reset_plot();
            for(int ptcl = 0; ptcl < path->getParam()->getNumParticles(); ptcl++){
                g1.plot_x(path->getBeads()->returnLinkedList(ptcl), " " );
            }*/
        }
    }
    
    std::cout << "\nCenter of mass acceptance: "<< 1.0*numacceptc/(numSteps) << "\nStaging acceptance: "<< 1.0*numaccepts/(numSteps)<< "\nBisection acceptance: "<< 1.0*numacceptb/(numSteps) << "\n\n";
    
    return energytr;
}