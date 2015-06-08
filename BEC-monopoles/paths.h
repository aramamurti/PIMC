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
#include "parameters.h"

using namespace std;

class paths{
public:
    paths(int procnum);
    ~paths();
    double vext(int slice, int ptcl);
    double potentialAction(int slice);
    double kineticAction(int slice, int dist);
    double kineticEnergy();
    double potentialEnergy();
    double energy();
    double cv();
    parameters* getParam(){return param;}
    void constPerms();
    void recompSingProb(vector<int> chdpart, int stslice);
    int factorial(int n){return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;}
    int permutation(int n, int k){return factorial(n)/factorial(n-k);}
    vector<vector<double>>* getProbList(){return &probList;}
    vector<vector<int>>* getPermList(){return &permList;}
    void setLast(vector<int> chdpart){last_chgd_part = chdpart;}
    vector<int> getLastP(){return last_chgd_part;}
    int getDist(){return multistep_dist;}
    vector<int> getNextConnection(){return nextConnection;}
    void setNextConnection(int ptcl, int newlink){nextConnection[ptcl] = newlink;}
    void setNextConnection(vector<int> nConn){nextConnection = nConn;}
    utility* getUte(){return ute;}
    
    void print();
    
    vector<vector<vector<double>>> beads;


    

private:
    parameters* param;
    potentials* pot;
    vector<vector<int>> permList;
    vector<vector<int>> permPart;
    vector<vector<double>> probList;
    utility* ute;
    vector<int> last_chgd_part;
    int multistep_dist;
    vector<int> nextConnection;
    bool printed;

};

#endif /* defined(__PIMCtest__paths__) */
