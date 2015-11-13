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
    void set_up_beads();
    void put_in_box();
    
    //getter methods
    int get_multistep_dist(){return multistep_dist;}
    boost::shared_ptr<Utility> get_util(){return util;}
    boost::shared_ptr<Parameters> get_parameters(){return params;}
    list_ptr get_beads(){return beads;}
    iVector get_last_changed(){return last_chd_parts;}
    iVector get_last_start_end(){iVector se; se.push_back(last_start); se.push_back(last_end); return se;}
    void set_last_changed(iVector lc){last_chd_parts = lc;}
    void set_last_start_end(int s, int e){last_start = s; last_end = e;}
    int get_processor_num(){return pnum;}
    iVector get_charge_list(){return charge_list;}
    
    float multvec[4];
    
    
private:
    
    class Separation_Table{
        
    public:
        
        Separation_Table(boost::shared_ptr<Path> path);
        ~Separation_Table(){};
        
        void set_up_st(){};
        void update_table();
        void add_bead();
        void remove_bead();
        
    private:
        
        boost::shared_ptr<Path> path;
        boost::unordered_map<int, std::pair<int, int>> bead_pairs;
        boost::unordered_map<int, std::vector<float> > bead_pair_sep;
        
        float rcut;
        
    };
    
    boost::shared_ptr<Parameters> params;
    boost::shared_ptr<Potential_Functions> pot;
    boost::shared_ptr<Utility> util;
    
    list_ptr beads;
    
    iVector charge_list;
    
    iVector last_chd_parts;
    int last_start;
    int last_end;
    
    int multistep_dist;
    int pnum;
    
    boost::shared_ptr<Separation_Table> sep_table;
    
};



#endif /* defined(__PIMCtest__paths__) */
