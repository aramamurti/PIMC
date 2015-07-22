//
//  potentials.cpp
//  PIMC
//
//  Created by Adith Ramamurti on 5/14/15.
//  Copyright (c) 2015 Adith Ramamurti. All rights reserved.
//

#include "potentials.h"

double potentials::harmonicPotential(std::vector<double> loc, double m, double w){
    double potVal = 0.0;
    for(int i = 0; i < loc.size(); i++){
        potVal += 0.5*m*pow(w,2)*pow(loc[i],2);
    }
    return potVal;
}

double potentials::lj_int(double dist){
    double eps = 10.9;
    double sig = 2.64;
    return 4*eps*(pow((sig/dist),12) - pow((sig/dist),6));
}

double potentials::hardSphere(double dist){
    if(dist < 140E-12)
        return 10E20;
    else
        return 0;
}

std::vector<double> potentials::grad_lj(std::vector<double> distVec, double dist){
    double eps = 10.9;
    double sig = 2.64;
    double dV = 24*eps*pow(sig,6)*(pow(dist,6)-2*pow(sig,6))/pow(dist,14);
    std::vector<double> gradV;
    for(std::vector<double>::iterator it = distVec.begin(); it != distVec.end(); it++){
        gradV.push_back(*it*dV);
    }
    
    return gradV;
}