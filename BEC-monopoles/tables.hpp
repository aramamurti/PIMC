//
//  tables.hpp
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


class Neighbor_Table{
    
public:
    
    Neighbor_Table(boost::shared_ptr<Path> path, float cutoff){};
    ~Neighbor_Table(){};
    
    void set_up_nn();
    void update_table();
    void add_bead();
    void remove_bead();

private:
    
    boost::shared_ptr<Path> path;
    boost::unordered_map<int, int> bead_grid;
    
    float rcut;

    
};


class Permutation_Table{
    
public:
    
    Permutation_Table(boost::shared_ptr<Path> path, float cutoff = -1);
    ~Permutation_Table(){};
    
    void set_up_perms();
    float recalc_perms(iVector ptcls, int slice);
    iVector pick_permutation(int ptcl, int start);
    
    iiVector* get_perm_list(){return &perm_list;}
    ffVector* get_prob_list(){return &prob_list;}
    
private:
    
    boost::shared_ptr<Path> path;
    
    bool cutoff;
    int multistep_dist;
    
    boost::shared_ptr<Neighbor_Table> ntable;
    
    iiVector perm_list;
    ffVector prob_list;
    iiVector perm_part_loc;
    iiVector permed_parts;

    
};


#endif /* tables_hpp */
