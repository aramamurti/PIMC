//
//  paths.h
//  PIMCtest
//
//  Created by Adith Ramamurti on 5/20/15.
//  Copyright (c) 2015 Adith Ramamurti. All rights reserved.
//

#ifndef __PIMCtest__paths__
#define __PIMCtest__paths__

#include <stdio.h>
#include <vector>
#include <set>
#include "potentials.h"
#include "utility.h"


using namespace std;

class paths{
public:
    paths(vector<vector<double>> beads, double tau, double lam, bool boson);
    ~paths();
    double vext(double R);
    double potentialAction(int slice);
    double kineticAction(int slice, int dist);
    double kineticEnergy();
    double potentialEnergy();
    double energy();
    int getNumParticles();
    int getNumSlices();
    double getTau();
    double getLam();
    void constPerms();
    int factorial(int n){return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;}
    int permutation(int n, int k){return factorial(n)/factorial(n-k);}
    vector<vector<int>> getPermList();
    utility* getUte();
    bool isBoson();
    
    vector<vector<double>> beads;

    

private:
    double tau;
    double lam;
    int numTimeSlices;
    int numParticles;
    potentials* pot;
    vector<vector<int>> permList;
    utility* ute;
    bool boson;

};

#endif /* defined(__PIMCtest__paths__) */
