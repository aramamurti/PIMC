////
////  EnergyEstimators.cpp
////  PIMC
////
////  Created by Adith Ramamurti on 5/13/15.
////  Copyright (c) 2015 Adith Ramamurti. All rights reserved.
////
//
//#include "energyEstimators.h"
//#include "potentials.h"
//#include <iostream>
//#include <vector>
//#include <cmath>
//#include "constants.h"
//#include "paths.h"
//
//using namespace std;
//
//energyEstimators::energyEstimators(){
//    unique_ptr<potentials> pot(new potentials());
//   
//    }
//energyEstimators::~energyEstimators(){
//    
//}
//
////double energyEstimators::thermodynamicEnergy(vector<vector<double>> beads){
////    double energy = 0.0;
////    int M = (int)beads.size();
////    double tau = constants().getTau();
////    double m = constants().getM();
////    double w = constants().getOmega();
////    double lambda = constants().getLambda();
////    int ndim = constants().getNdim();
////    
////    double t1 = ndim*0.5/tau;
////    double t2 = 0.0, t3 = 0.0;
////    
////    for(int i = 0; i < M; i++){
////        for(int j = 0; j < ndim; j++){
////            t2 -= pow((beads[i][j]-beads[(i+1)%M][j]),2);
////        }
////        
////        t3 += tau/2.0*(pot->harmonicPotential(beads[i],m,w)-pot->harmonicPotential(beads[(i+1)%M],m,w));
////    }
////    t2 *= 1/(4*lambda*pow(tau,2)*M);
////    t3 /= M;
////    
////    energy = t1 + t2 + t3;
////    return energy;
////}
//
//
//
////double energyEstimators::virialEnergy(vector<vector<double>> beads){
////    double energy = 0.0;
////    int M = (int) beads.size();
////    int virialWindow = M;
////    double tau = constants().getTau();
////    double m = constants().getM();
////    double w = constants().getOmega();
////    double lambda = constants().getLambda();
////    int ndim = constants().getNdim();
////    
////    double t1 = ndim*0.5/(tau*virialWindow);
////    
////    double t2 = 0.0, t3 = 0.0, t4 = 0.0;
////    for(int i = 0; i < M; i++){
////        
////        for(int j = 0; j < ndim; j++){
////            t2 -= (beads[i][j]-beads[(i+virialWindow)%M][j])*(beads[i][j]-beads[(i+1)%M][j]);
////        }
////        
////        double spadelt = 0.0;
////        for(int k = -virialWindow+1; k < virialWindow; k++){
////            for(int j = 0; j < ndim; j++){
////                spadelt += beads[i][j] - beads[(i+M+k)%M][j];
////            }
////        }
////
////        t3 += (pot->harmonicPotentialDeriv(beads[(i-1+M)%M], m, w)+pot->harmonicPotentialDeriv(beads[i], m,w))*spadelt;
////        t4 += pot->harmonicPotential(beads[i], m, w);
////    }
////    t2 /= (4*virialWindow*lambda*pow(tau,2));
////    t3 /= 4.0*virialWindow;
////    energy = t1 + (t2 + t3 + t4)/M;
////    return energy;
////    
////}
//
//
