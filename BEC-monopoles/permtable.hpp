//
//  permtable.hpp
//  BEC-monopoles
//
//  Created by Adith Ramamurti on 11/7/15.
//  Copyright © 2015 Adith Ramamurti. All rights reserved.
//

#ifndef tables_hpp
#define tables_hpp

#include <stdio.h>
#include "uni_header.h"
#include "path.h"
#include "actions.hpp"

class Permutation_Table{
    
public:
    
    Permutation_Table(boost::shared_ptr<Path> path);
    ~Permutation_Table(){};
    
    void set_up_perms();
    double recalc_perms(iVector ptcls, int slice);
    iVector pick_permutation(int ptcl, int start);
    
    iiVector* get_perm_list(){return &perm_list;}
    ddVector* get_prob_list(){return &prob_list;}
    
private:
    
    boost::shared_ptr<Path> path;
    
    int multistep_dist;
    
    iiVector perm_list;
    ddVector prob_list;
    iiVector perm_part_loc;
    iiVector permed_parts;
    
    boost::shared_ptr<Kinetic_Action> ka;

    
};


#endif /* tables_hpp */
