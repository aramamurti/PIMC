//
//  main.cpp
//  PIMC
//
//  Created by Adith Ramamurti on 5/20/15.
//  Copyright (c) 2015 Adith Ramamurti. All rights reserved.
//


/***********************/
 
#include <iostream>
#include <vector>
#include "moves.h"
#include "paths.h"
#include "pimc.h"


using namespace std;


int main(int argc, const char * argv[]) {
    double T = 0.05;
    double lam = 0.5;
    
    int numParticles = 1;
    int numTimeSlices = 40;
    int numSteps = 300000;
    double tau = 1/(T*numTimeSlices);
    
    cout<< "Simulation Parameters:\n";
    cout<< "N      = \t" << numParticles <<"\n";
    cout<< "tau    = \t" << tau << "\n";
    cout << "lambda =\t" << lam <<"\n";
    cout<< "T      = \t" << T << "\n\n";
        
    vector<vector<double>> beads(numTimeSlices, vector<double>(numParticles,0.0));
    
    paths* path = new paths(beads, tau, lam, true);
    pimc sim;
    vector<double> energy = sim.run(numSteps, path);
    delete path;
    
    cout<< "Energy = " <<utility().vecavg(energy) << " +/- "<< utility().vecstd(energy)/sqrt(energy.size())<<"\n";
    
    return 0;
}

