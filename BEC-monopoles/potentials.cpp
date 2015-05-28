//
//  potentials.cpp
//  PIMC
//
//  Created by Adith Ramamurti on 5/14/15.
//  Copyright (c) 2015 Adith Ramamurti. All rights reserved.
//

#include "potentials.h"
#include <cmath>
#include <vector>
#include <iostream>
#include "constants.h"


using namespace std;

double potentials::harmonicPotential(double pos, double m, double w){
    double potVal = 0.0;
    potVal += 0.5*m*pow(w,2)*pow(pos,2);
    return potVal;
}

double potentials::harmonicPotentialDeriv(vector<double> pos, double m, double w){
    int ndim = constants().getNdim();
    double hpdiv = 0;
    for(int j = 0; j < ndim; j++){
        hpdiv += pos[j];
    }
    hpdiv *= m*w;
    return hpdiv;
}

double potentials::harmonicPotentialVir(vector<double> pos,double m, double w){
    int ndim = constants().getNdim();
    double potVal = 0.0;
    for(int j =0; j< ndim; j++) {
        potVal += pos[j]*m*pow(w,2)*pos[j];
    }
    return potVal;
}