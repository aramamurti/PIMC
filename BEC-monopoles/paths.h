//
//  paths.h
//  PIMCtest
//
//  Created by Adith Ramamurti on 5/20/15.
//  Copyright (c) 2015 Adith Ramamurti. All rights reserved.
//

#ifndef __PIMCtest__paths__
#define __PIMCtest__paths__

#include "potentials.h"
#include "utility.h"
#include "parameters.h"
#include "container.h"

class paths{
public:
    typedef std::tr1::shared_ptr<LinkedList<std::vector<double>>> list_ptr;

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
    double recompSingProb(int stslice);
    
    //getter methods
    int getDist(){return multistep_dist;}
    utility* getUte(){return ute;}
    parameters* getParam(){return param;}
    std::vector<std::vector<double>>* getProbList(){return &probList;}
    std::vector<std::vector<int>>* getPermList(){return &permList;}
    list_ptr getBeads(){return beads;}

    
    bool printed;

private:
    

    parameters* param;
    potentials* pot;
    list_ptr beads;
    std::vector<std::vector<int>> permList;
    std::vector<std::vector<int>> permPart;
    std::vector<std::vector<double>> probList;
    utility* ute;
    int multistep_dist;
    double multvec[4];


};

#endif /* defined(__PIMCtest__paths__) */
