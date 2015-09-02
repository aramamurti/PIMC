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
    double eps = 10.22;
    double sig = 2.556;
    return 4*eps*(pow((sig/dist),12) - pow((sig/dist),6));
}

double potentials::hardSphere(double dist){
    if(dist < 140E-12)
        return 10E20;
    else
        return 0;
}

double potentials::aziz_int(double dist){
    double val = 10.8*(544850.4 * exp(-4.50018*dist)-(9424.94/pow(dist,10)+2556.63/pow(dist,8)+937.38/pow(dist,6)) *aziz_pcws(dist));
    return val;
}

inline double potentials::aziz_pcws(double dist){
    if(dist >= 3.68335)
        return 1;
    else
        return exp(-pow((3.68335/dist-1),2));
}