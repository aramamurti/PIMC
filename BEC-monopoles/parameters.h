//
//  parameters.h
//  BEC-monopoles
//
//  Created by Adith Ramamurti on 5/29/15.
//  Copyright (c) 2015 Adith Ramamurti. All rights reserved.
//

#ifndef __BEC_monopoles__parameters__
#define __BEC_monopoles__parameters__

#include <stdio.h>
#include <cmath>

using namespace std;

class parameters{
    
private:
    double tau;
    double T;
    double hbar;
    double kb;
    double w;
    double m;
    double lambda;
    int ndim;
    double lam;
    int numTimeSlices;
    int numParticles;
    bool boson;
    int numSteps;
    int skip;
    int equil;
    double boxsize;
    bool pbc;
    vector<bool> pmv;
    
public:
    parameters(){
        ndim = 1;
        m = 1, hbar = 1.0, kb = 1.0, w = 1.0;
        T = 0.2;
        lam = pow(hbar,2)/(2*m);
        
        boson = true;
        
        numParticles = 3;
        numTimeSlices = 80;
        numSteps = 300000;
        skip = 100;
        equil = 10000;
        
        tau = 1/(T*numTimeSlices);
        pmv = {true, true, false};
        
        pbc = false;
        if(pbc)
            boxsize = 6e-10;
        else
            boxsize = -1;

    }
    ~parameters(){}
    
    void setT(double newT){
        T = newT;
        tau = 1/(T*numTimeSlices);
    }
    int getndim(){return ndim;}
    double getT(){return T;}
    double gethbar(){return hbar;}
    double getkb(){return kb;}
    double getomega(){return w;}
    double getm(){return m;}
    bool getboson(){return boson;}
    double gettau(){return tau;}
    double getlam(){return lam;}
    int getNumSteps(){return numSteps;}
    int getSkip(){return skip;}
    int getEquil(){return equil;}
    int getNumTimeSlices(){return numTimeSlices;}
    int getNumParticles(){return numParticles;}
    double getBoxSize(){return boxsize;}
    vector<bool> getMoves(){return pmv;}
};

#endif /* defined(__BEC_monopoles__parameters__) */
