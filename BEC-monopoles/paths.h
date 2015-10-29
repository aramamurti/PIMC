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
    paths(int procnum, std::ofstream &f);
    ~paths();
    
    //methods
    double vext(int slice, int ptcl);
    double potentialAction(int slice);
    double kineticAction(int slice, int dist);
    double relativisticKineticAction(int slice, int dist);
    double kineticEnergy();
    double relativisticKineticEnergy();
    double potentialEnergy();
    double energy();
    std::vector<int> getCycles();
    std::vector<int> getWindingNumber();
    double virialEnergy();
    double relenergy();
    double cv();
    void constPerms(int procnum);
    double recompSingProb(std::vector<int> ptcls, int stslice);
    void putInBox();
    
    //getter methods
    int getDist(){return multistep_dist;}
    utility* getUte(){return ute;}
    parameters* getParam(){return param;}
    std::vector<std::vector<double>>* getProbList(){return &probList;}
    std::vector<std::vector<int>>* getPermList(){return &permList;}
    list_ptr getBeads(){return beads;}
    std::vector<int> getLastChgd(){return lastChdParticles;}
    std::vector<int> getlastLocs(){std::vector<int> locs; locs.push_back(laststart); locs.push_back(lastend); return locs;}
    void setlastChgd(std::vector<int> lc){lastChdParticles = lc;}
    void slstep(int s, int e){laststart = s; lastend = e;}
    
    int getPNum(){return pnum;}
    
    
    int numswap;

    
private:
    

    parameters* param;
    potentials* pot;
    list_ptr beads;
    std::vector<std::vector<int>> permList;
    std::vector<std::vector<int>> permPart;
    std::vector<std::vector<int>> permPartLoc;
    std::vector<std::vector<double>> probList;
    
    std::vector<int> lastChdParticles;
    int laststart;
    int lastend;
    
    utility* ute;
    int multistep_dist;
    double multvec[4];
    int pnum;


};

#endif /* defined(__PIMCtest__paths__) */
