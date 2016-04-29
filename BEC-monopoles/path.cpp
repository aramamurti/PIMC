//
//  paths.cpp
//  PIMCtest
//
//  Created by Adith Ramamurti on 5/20/15.
//  Copyright (c) 2015 Adith Ramamurti. All rights reserved.
//

#include "path.h"

/****
 Constructor: sets up paths and parameters for each particle and constructs all objects necessary for the PIMC simulation steps.
 *****/
Path::Path(int procnum, IO &writer, boost::shared_ptr<Parameters> parameters)
:multvec{1.0,20.0,100.0, 400.0}{ // Multiplication factor for the permute probabilities
    
    //Declare all objects and variables, set all parameters that need to be modified from default (in parameters.h)
    
    util = boost::shared_ptr<Utility>(new Utility(procnum));
    params = parameters;
    
    writer.write_parameters(params);
    multistep_dist = 8;
    last_start = 0;
    last_end = 0;
    pnum = procnum;
    
    set_up_beads();
    
}

//Deconstructor
Path::~Path(){
    beads.reset();
}

void Path::set_up_beads(){
    //Set up the initial distribution of particles.
    ddVector offset(params->get_num_particles(), dVector(params->get_ndim(), 0.0));
    
    if(params->get_box_size() == -1){
        for(unsigned int ptcl = 0; ptcl < params->get_num_particles(); ptcl++){
            for(unsigned int ndim = 0; ndim < params->get_ndim(); ndim++){
                offset[ptcl][ndim] = util->randnormed(1)-0.5;
            }
        }
        
    }
    else{
        offset.resize(0);
        for(int i = 0; i < params->get_num_particles(); i++){
            dVector pos(0);
            for(int ndim = 0; ndim < params->get_ndim(); ndim++){
                pos.push_back(util->randnormed(params->get_box_size()));
            }
            offset.push_back(pos);
        }
    }
    
    beads = list_ptr(new PathList<dVector>(params->get_ndim(), multistep_dist, params->get_box_size()));
    
    for(int slice = 0; slice < params->get_num_timeslices(); slice++){
        for(int ptcl = 0; ptcl < params->get_num_particles(); ptcl++){
            beads->push_back(offset[ptcl],ptcl);
        }
    }
    
    beads->generate_neighbors();
    beads->generate_perm_seps();
    
    
    if(params->is_charged()){
        int numcharges = params->get_num_chgs();
        for(int i = 0; i < params->get_num_particles();i++){
            charge_list.push_back(2*(i%numcharges)-1);
        }
    }
    
    for(int i = 0; i < charge_list.size(); i++){
        beads->set_charge(i, charge_list[i]);
    }
    
    //Make the beads a periodic chain
    beads->make_circular();
}

void Path::put_in_box(){
    
    if(params->get_box_size() != -1)
        for(int ptcl = 0; ptcl < beads->get_num_particles(); ptcl++)
            for(int slice = 0; slice < params->get_num_timeslices(); slice++){
                beads->set_bead_data(ptcl, slice, util->location(beads->get_bead_data(ptcl, slice), params->get_box_size()),beads->get_bead_data(ptcl, slice));
            }
    
}

void Path::put_in_box(iVector changed_ptcls, int start_slice, int end_slice){
    
    if(params->get_box_size() != -1)
        for(std::vector<int>::iterator it = changed_ptcls.begin(); it != changed_ptcls.end(); it++)
            for(int slice = start_slice; slice < end_slice; slice++){
                beads->set_bead_data(*it, slice, util->location(beads->get_bead_data(*it, slice), params->get_box_size()),beads->get_bead_data(*it, slice));
            }
    
}