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


paths::paths(int procnum){
    ute = new utility(procnum);
    param = new parameters();
    
    cout<< "Simulation Parameters:\n" << "N      = \t" << param->getNumParticles() <<"\n" <<"tau    = \t" << param->gettau() << "\n" << "lambda =\t" << param->getlam() <<"\n" << "T      = \t" << param->getT() << "\n\n";
    
    vector<vector<vector<double>>> beads(param->getNumTimeSlices(), vector<vector<double>>(param->getNumParticles(),vector<double>(param->getndim(), 0.0)));
    this->beads = beads;
    for(int slice = 0; slice < param->getNumTimeSlices(); slice++){
        for(int ptcl = 0; ptcl < param->getNumParticles(); ptcl++){
            for(int ndim = 0; ndim < param->getndim(); ndim++){
                this->beads[slice][ptcl][ndim]= 0.5*(ute->randnormed(2)-1);
            }
        }
    }

    pot = new potentials();
    constPerms();
}

paths::~paths(){
    delete pot;
    delete ute;
    delete param;
}


void paths::constPerms(){
    vector<int> d(param->getNumParticles());
    int k;
    if(param->getNumParticles() <= 4)
        k = param->getNumParticles();
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
        vector<int> identity(param->getNumParticles());
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
    for(int ptcl = 0; ptcl < param->getNumParticles(); ptcl++){
        for(int ndim = 0; ndim < param->getndim(); ndim++){
            pot += vext(beads[slice][ptcl][ndim]);
        }
    }
    return param->gettau()*pot;
}

double paths::kineticAction(int slice, int dist){
    double kin = 0;
    for(int ptcl = 0; ptcl < param->getNumParticles(); ptcl++){
        for(int ndim = 0; ndim < param->getndim(); ndim++){
            kin += 1/(2*param->getNumTimeSlices()*dist*param->gettau())*pow(beads[slice][ptcl][ndim]-beads[(slice-dist+param->getNumTimeSlices())%param->getNumTimeSlices()][ptcl][ndim],2);
        }
    }
    return kin;
}

double paths::potentialEnergy(){
    double PE = 0.0;
    double m = param->getm();
    double w = param->getomega();
    for(int slice = 0; slice<param->getNumTimeSlices();slice++){
        for(int ptcl = 0; ptcl<param->getNumParticles(); ptcl++){
            for(int ndim = 0; ndim < param->getndim(); ndim++){
                double R = beads[slice][ptcl][ndim];
                PE += pot->harmonicPotential(R, m, w);
            }
        }
    }
    PE = PE/param->getNumTimeSlices();
    return PE;
}

double paths::kineticEnergy(){
    double tot = 0.0;
    double norm = 1.0/(4.0*param->getlam()*pow(param->gettau(),2));
    for(int slice = 0; slice < param->getNumTimeSlices(); slice++){
        int slicep1 = (slice+1)%param->getNumTimeSlices();
        for(int ptcl = 0; ptcl < param->getNumParticles(); ptcl++){
            for(int ndim = 0; ndim < param->getndim(); ndim++){
                double delR = beads[slicep1][ptcl][ndim]-beads[slice][ptcl][ndim];
                tot -= norm*delR*delR;
            }
        }
    }
    double KE = 0.5*param->getndim()*param->getNumParticles()/param->gettau() +tot/param->getNumTimeSlices();
    return KE;
}


double paths::energy(){
    double energy = kineticEnergy()+potentialEnergy();
    return energy;
}

vector<vector<int>> paths::getPermList(){
    return permList;
}

utility* paths::getUte(){
    return ute;
}

parameters* paths::getParam(){
    return param;
}