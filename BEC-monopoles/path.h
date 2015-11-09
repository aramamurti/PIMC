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
#include "IO.hpp"


class Path{
public:
    typedef boost::shared_ptr<PathList<vectorf>> list_ptr;

    //constructor and destructor
    Path(int procnum, IO &writer);
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
    utility* getUte(){return util;}
    Parameters* get_parameters(){return params;}
    vectorff* get_prob_list(){return &prob_list;}
    vectorii* get_perm_list(){return &perm_list;}
    list_ptr get_beads(){return beads;}
    vectori get_last_changed(){return last_chd_parts;}
    vectori get_last_locs(){std::vector<int> locs; locs.push_back(last_start); locs.push_back(last_end); return locs;}
    void set_last_changed(vectori lc){last_chd_parts = lc;}
    void set_last_step(int s, int e){last_start = s; last_end = e;}
    int getPNum(){return pnum;}
    
private:
    

    Parameters* params;
    potentials* pot;
    list_ptr beads;
    vectorii perm_list;
    vectorii permed_parts;
    vectorii perm_part_loc;
    vectorff prob_list;
    
    vectori charge_list;
    
    vectori last_chd_parts;
    int last_start;
    int last_end;
    
    utility* util;
    int multistep_dist;
    float multvec[4];
    int pnum;


};

#endif /* defined(__PIMCtest__paths__) */
