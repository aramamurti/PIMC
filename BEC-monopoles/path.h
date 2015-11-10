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
    typedef boost::shared_ptr<PathList<fVector> > list_ptr;

    //constructor and destructor
    Path(int procnum, IO &writer);
    ~Path();
    
    //methods
    float vext(int slice, int ptcl);
    float potential_action(int slice);
    float kinetic_action(int slice, int dist);
    float kinetic_energy();
    float potential_energy();
    float energy();
    iVector get_cycles();
    iVector get_winding_number();
    void set_up_beads();
    void constr_perms(int procnum);
    float slice_perm_prob(iVector ptcls, int stslice);
    void put_in_box();
    
    //getter methods
    int get_multistep_dist(){return multistep_dist;}
    boost::shared_ptr<Utility> get_util(){return util;}
    boost::shared_ptr<Parameters> get_parameters(){return params;}
    ffVector* get_prob_list(){return &prob_list;}
    iiVector* get_perm_list(){return &perm_list;}
    list_ptr get_beads(){return beads;}
    iVector get_last_changed(){return last_chd_parts;}
    iVector get_last_locs(){iVector locs; locs.push_back(last_start); locs.push_back(last_end); return locs;}
    void set_last_changed(iVector lc){last_chd_parts = lc;}
    void set_last_step(int s, int e){last_start = s; last_end = e;}
    int getPNum(){return pnum;}
    
private:
    

    boost::shared_ptr<Parameters> params;
    boost::shared_ptr<potentials> pot;
    boost::shared_ptr<Utility> util;

    list_ptr beads;
    iiVector perm_list;
    iiVector permed_parts;
    iiVector perm_part_loc;
    ffVector prob_list;
    
    iVector charge_list;
    
    iVector last_chd_parts;
    int last_start;
    int last_end;
    
    int multistep_dist;
    float multvec[4];
    int pnum;


};

#endif /* defined(__PIMCtest__paths__) */
