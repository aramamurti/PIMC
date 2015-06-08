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

class paths{
public:
    
    //constructor and destructor
    paths(int procnum);
    ~paths();
    
    //methods
    double vext(int slice, int ptcl);
    double potentialAction(int slice);
    double kineticAction(int slice, int dist);
    double kineticEnergy();
    double potentialEnergy();
    double energy();
    double cv();
    void constPerms();
    void recompSingProb(std::vector<int> chdpart, int stslice);
    
    //setter methods
    void setLast(std::vector<int> chdpart){last_chgd_part = chdpart;}
    void setNextConnection(int ptcl, int newlink){nextConnection[ptcl] = newlink;}
    void setNextConnection(std::vector<int> nConn){nextConnection = nConn;}
    
    
    //getter methods
    int getDist(){return multistep_dist;}
    utility* getUte(){return ute;}
    parameters* getParam(){return param;}
    std::vector<int> getLastP(){return last_chgd_part;}
    std::vector<int> getNextConnection(){return nextConnection;}
    std::vector<std::vector<double>>* getProbList(){return &probList;}
    std::vector<std::vector<int>>* getPermList(){return &permList;}

    void print();
    
    std::vector<std::vector<std::vector<double>>> beads;


    

private:
    parameters* param;
    potentials* pot;
    std::vector<std::vector<int>> permList;
    std::vector<std::vector<int>> permPart;
    std::vector<std::vector<double>> probList;
    utility* ute;
    std::vector<int> last_chgd_part;
    int multistep_dist;
    std::vector<int> nextConnection;
    bool printed;

};

#endif /* defined(__PIMCtest__paths__) */
