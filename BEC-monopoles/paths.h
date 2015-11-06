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
    typedef std::tr1::shared_ptr<LinkedList<vectorf>> list_ptr;

    //constructor and destructor
    paths(int procnum, std::ofstream &f);
    ~paths();
    
    //methods
    float vext(int slice, int ptcl);
    float potentialAction(int slice);
    float kineticAction(int slice, int dist);
    float kineticEnergy();
    float potentialEnergy();
    float energy();
    vectori getCycles();
    vectori getWindingNumber();
    void setup_beads();
    void constr_perms(int procnum);
    float slice_perm_prob(vectori ptcls, int stslice);
    void putInBox();
    
    //getter methods
    int getDist(){return multistep_dist;}
    utility* getUte(){return ute;}
    parameters* getParam(){return param;}
    vectorff* getProbList(){return &probList;}
    vectorii* getPermList(){return &permList;}
    list_ptr getBeads(){return beads;}
    vectori getLastChgd(){return lastChdParticles;}
    vectori getlastLocs(){std::vector<int> locs; locs.push_back(last_start); locs.push_back(last_end); return locs;}
    void setlastChgd(vectori lc){lastChdParticles = lc;}
    void slstep(int s, int e){last_start = s; last_end = e;}
    int getPNum(){return pnum;}
    
private:
    

    parameters* param;
    potentials* pot;
    list_ptr beads;
    vectorii permList;
    vectorii permPart;
    vectorii permPartLoc;
    vectorff probList;
    
    vectori chgList;
    
    vectori lastChdParticles;
    int last_start;
    int last_end;
    
    utility* ute;
    int multistep_dist;
    float multvec[4];
    int pnum;


};

#endif /* defined(__PIMCtest__paths__) */
