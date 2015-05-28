//
//  paths.cpp
//  PIMCtest
//
//  Created by Adith Ramamurti on 5/20/15.
//  Copyright (c) 2015 Adith Ramamurti. All rights reserved.
//

#include "paths.h"
#include <iostream>
#include <cmath>
#include <numeric>
#include <algorithm>
#include "constants.h"


paths::paths(vector<vector<double>> beads, double tau, double lam, bool boson){
    ute = new utility();
    
    this->boson = boson;
    this->tau = tau;
    this->lam = lam;
    this->beads = beads;
    numTimeSlices = (int) beads.size();
    numParticles = (int) beads[0].size();
    for(int slice = 0; slice < numTimeSlices; slice++){
        for(int ptcl = 0; ptcl < numParticles; ptcl++){
            this->beads[slice][ptcl]= 0.5*(ute->randnormed(2)-1);
        }
    }
    unique_ptr<potentials> pot(new potentials());
    constPerms();
}

paths::~paths(){
    delete pot;
    delete ute;
}


void paths::constPerms(){
    vector<int> d(numParticles);
    int k;
    if(numParticles <= 4)
        k = numParticles;
    else
        k = 4;
    
    vector<vector<int>> initPermList = {};
    iota(d.begin(),d.end(),0);

    do
    {
        vector<int> tempvec = {};
        for (int i = 0; i < k; i++)
        {
            tempvec.push_back(d[i]);
        }
        initPermList.push_back(tempvec);
        std::reverse(d.begin()+k,d.end());
    } while (next_permutation(d.begin(),d.end()));
    
    for(int i = 0; i < initPermList.size(); i++){
        vector<int> identity(numParticles);
        iota(identity.begin(),identity.end(),0);
        int displaced[k];
        for(int j = 0; j < k; j++){
            displaced[j] = initPermList[i][j];
        }
        sort(initPermList[i].begin(),initPermList[i].end());
        for(int j = 0; j < k; j++){
            identity[initPermList[i][j]] = displaced[j];
        }
        permList.push_back(identity);
    }
    
    sort(permList.begin(), permList.end());
    auto last = unique(permList.begin(), permList.end());
    permList.erase(last, permList.end());
}

double paths::vext(double R){
    return 0.5*R*R;
}

double paths::potentialAction(int slice){
    double pot = 0;
    for(int ptcl = 0; ptcl < numParticles; ptcl++){
        pot += vext(beads[slice][ptcl]);
    }
    return tau*pot;
}

double paths::kineticAction(int slice, int dist){
    double kin = 0;
    for(int ptcl = 0; ptcl < numParticles; ptcl++){
        kin += 1/(2*numTimeSlices*dist*tau)*pow(beads[slice][ptcl]-beads[(slice-dist+numTimeSlices)%numTimeSlices][ptcl],2);
    }
    return kin;
}

double paths::potentialEnergy(){
    double PE = 0.0;
    double m = constants().getM();
    double w = constants().getOmega();
    for(int slice = 0; slice<getNumSlices();slice++){
        for(int ptcl = 0; ptcl<getNumParticles(); ptcl++){
            double R = beads[slice][ptcl];
            PE += pot->harmonicPotential(R, m, w);
        }
    }
    PE = PE/getNumSlices();
    return PE;
}

double paths::kineticEnergy(){
    double tot = 0.0;
    double norm = 1.0/(4.0*getLam()*pow(getTau(),2));
    for(int slice = 0; slice < getNumSlices(); slice++){
        int slicep1 = (slice+1)%getNumSlices();
        for(int ptcl = 0; ptcl < getNumParticles(); ptcl++){
            double delR = beads[slicep1][ptcl]-beads[slice][ptcl];
            tot -= norm*delR*delR;
        }
    }
    double KE = 0.5*getNumParticles()/getTau() +tot/getNumSlices();
    return KE;
}

double paths::energy(){
    double energy = kineticEnergy()+potentialEnergy();
    return energy;
}

int paths::getNumParticles(){
    return numParticles;
}
int paths::getNumSlices(){
    return numTimeSlices;
}

double paths::getTau(){
    return tau;
}

double paths::getLam(){
    return lam;
}

vector<vector<int>> paths::getPermList(){
    return permList;
}

utility* paths::getUte(){
    return ute;
}

bool paths::isBoson(){
    return boson;
}