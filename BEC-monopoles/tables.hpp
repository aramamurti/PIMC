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

class TableBase{
public:
    TableBase();
    ~TableBase();
    void set_up_table();
    
private:
    
};

class PermTable: public TableBase{
public:
    PermTable();
    ~PermTable();
    void set_up_perms();
    void recalc_perms();
    iiVector get_perm_list();
    ffVector get_perm_prob();
    
private:
    
    
    
};

class NNTable: public TableBase{
public:
    NNTable();
    ~NNTable();
    void set_up_nn();
    void update_table();
    void add_bead();
    void rem_bead();
    
private:
    
};

#endif /* tables_hpp */
