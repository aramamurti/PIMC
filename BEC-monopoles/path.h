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

class Path{
public:
    typedef std::tr1::shared_ptr<LinkedList<vectorf>> list_ptr;

    //constructor and destructor
    Path(int procnum, std::ofstream &f);
    ~Path();
    
    //methods
    float vext(int slice, int ptcl);
    float potentialAction(int slice);
    float kineticAction(int slice, int dist);
    float kineticEnergy();
    float potentialEnergy();
    float energy();
    vectori get_cycles();
    vectori get_winding_number();
    void setup_beads();
    void constr_perms(int procnum);
    float slice_perm_prob(vectori ptcls, int stslice);
    void put_in_box();
    
    //get_Ter methods
    int getDist(){return multistep_dist;}
    utility* getUte(){return ute;}
    Parameters* get_parameters(){return params;}
    vectorff* get_prob_list(){return &probList;}
    vectorii* get_perm_list(){return &permList;}
    list_ptr get_beads(){return beads;}
    vectori get_last_changed(){return lastChdParticles;}
    vectori get_last_locs(){std::vector<int> locs; locs.push_back(last_start); locs.push_back(last_end); return locs;}
    void set_last_changed(vectori lc){lastChdParticles = lc;}
    void set_last_step(int s, int e){last_start = s; last_end = e;}
    int getPNum(){return pnum;}
    
private:
    

    Parameters* params;
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
