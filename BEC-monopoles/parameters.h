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
    float tau;
    float T;
    float kb;
    float lambda;
    int ndim;
    float lam;
    int numTimeSlices;
    int numParticles;
    bool boson;
    int numSteps;
    int skip;
    int equil;
    float boxsize;
    bool pbc;
    bool pots[5];
    int numcharges;
    bool charged;

    
public:
    parameters(){
        ndim = 1;
        kb = 1.0;
        T = 1.0;
        lam = 0.5;//pow(hbar,2)/(2*m);
        
        numcharges = 2;
        charged = true;
        
        boson = false;
        
        numParticles = 6;
        numTimeSlices = 10;
        numSteps = 10000;
        skip = 10;
        equil = 1000;
        
        pots[0] = true;
        pots[1] = false;
        pots[2] = false;
        pots[3] = false;
        pots[4] = false;
        
        tau = 1/(T*numTimeSlices);
        
        pbc = false;
        if(pbc)
            boxsize = pow(numParticles/10.,1/3.)*7.7099;
        else
            boxsize = -1;

    }
    ~parameters(){}
    
    void setT(float newT){
        T = newT;
        tau = 1/(T*numTimeSlices);
    }
    void setNumTS(float newTS){
        numTimeSlices = newTS;
        tau = 1/(T*numTimeSlices);
    }
    
    int getndim(){return ndim;}
    float getT(){return T;}
    float getkb(){return kb;}
    bool isboson(){return boson;}
    bool ischarged(){return charged;}
    int getNumCharges(){return numcharges;}
    float gettau(){return tau;}
    float getlam(){return lam;}
    int getNumSteps(){return numSteps;}
    int getSkip(){return skip;}
    int getEquil(){return equil;}
    int getNumTimeSlices(){return numTimeSlices;}
    int getNumParticles(){return numParticles;}
    float getBoxSize(){return boxsize;}
    bool* getPots(){return pots;}
};

#endif /* defined(__BEC_monopoles__parameters__) */
