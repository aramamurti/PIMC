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
    vector<bool> pmv;
    
public:
    parameters(){
        bool input = false;
        
        ndim = 1;
        m = 1.0, hbar = 1.0, kb = 1.0, w = 1.0;
        T = 0.05;
        lam = 0.5;
        
        boson = true;
        
        numParticles = 1;
        numTimeSlices = 40;
        numSteps = 300000;
        tau = 1/(T*numTimeSlices);
        pmv = {true, true, false};
        
//        if(input){
//            cout << "Number of particles: ";
//            cin >> numParticles;
//            cout <<endl;
//            
//            cout << "Number of imaginary-time slices: ";
//            cin >> numTimeSlices;
//            cout <<endl;
//            
//            cout << "Are the particles bosons? (y/n::1/0) ";
//            cin >> boson;
//            cout <<endl;
//            
//            cout << "Number of imaginary-time slices: ";
//            cin >> numTimeSlices;
//            cout <<endl;
//            
//            cout << "Number of dimensions: ";
//            cin >> ndim;
//            cout << endl;
//            
//        }
        
        
    }
    ~parameters(){}
    
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
    int getNumTimeSlices(){return numTimeSlices;}
    int getNumParticles(){return numParticles;}
    vector<bool> getMoves(){return pmv;}
};

#endif /* defined(__BEC_monopoles__parameters__) */
