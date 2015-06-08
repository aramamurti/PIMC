//
//  potentials.cpp
//  PIMC
//
//  Created by Adith Ramamurti on 5/14/15.
//  Copyright (c) 2015 Adith Ramamurti. All rights reserved.
//

#include "potentials.h"
#include "parameters.h"
#include <cmath>
#include <vector>
#include <iostream>

double potentials::harmonicPotential(std::vector<double> loc, double m, double w){
    double potVal = 0.0;
    for(int i = 0; i < loc.size(); i++){
        potVal += 0.5*m*pow(w,2)*pow(loc[i],2);
    }
    return potVal;
}

double potentials::lj_int(double dist){
    double eps = 10.2;
    double sig = 2.28e-10;
    double potVal = 4*eps*(pow((sig/dist),12) - pow((sig/dist),6));
    return potVal;
}
