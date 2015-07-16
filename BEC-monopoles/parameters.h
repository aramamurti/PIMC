//
//  parameters.h
//  BEC-monopoles
//
//  Created by Adith Ramamurti on 5/29/15.
//  Copyright (c) 2015 Adith Ramamurti. All rights reserved.
//

#ifndef __BEC_monopoles__parameters__
#define __BEC_monopoles__parameters__

#include "uni_header.h"

class parameters{
    
private:
    double tau;
    double T;
    double kb;
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
    bool pots[2];

    
public:
    parameters(){
        ndim = 1;
        kb = 1.0;
        T = 0.2;
        lam = 0.5;//pow(hbar,2)/(2*m);
        
        boson = true;
        
        numParticles = 2;
        numTimeSlices = 80;
        numSteps = 100000;
        skip = 100;
        equil = 10000;
        pots[0] = false;
        pots[1] = true;
        
        tau = 1/(T*numTimeSlices);
        
        pbc = true;
        if(pbc)
            boxsize = 9.7;
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
    double getkb(){return kb;}
    bool getboson(){return boson;}
    double gettau(){return tau;}
    double getlam(){return lam;}
    int getNumSteps(){return numSteps;}
    int getSkip(){return skip;}
    int getEquil(){return equil;}
    int getNumTimeSlices(){return numTimeSlices;}
    int getNumParticles(){return numParticles;}
    double getBoxSize(){return boxsize;}
    bool* getPots(){return pots;}
};

#endif /* defined(__BEC_monopoles__parameters__) */
