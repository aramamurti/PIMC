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
    multistep_dist = 16;
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
//        offset.resize(0);
//        unsigned int pps = (int)ceil(pow(params->get_num_particles(),1.0/((double)params->get_ndim())));
//        double spacing = params->get_box_size()/pps;
//        for(int i = 0; i < pps; i++){
//            if(params->get_ndim() == 1){
//                dVector pos;
//                pos.push_back(i*spacing);
//                offset.push_back(pos);
//            }
//            else{
//                for(int j = 0; j < pps; j++){
//                    if(params->get_ndim() == 1){
//                        dVector pos;
//                        pos.push_back(i*spacing);
//                        pos.push_back(j*spacing);
//                        offset.push_back(pos);
//                    }
//                    else{
//                        for(int k = 0; k < pps; k++){
//                            dVector pos;
//                            pos.push_back(i*spacing);
//                            pos.push_back(j*spacing);
//                            pos.push_back(k*spacing);
//                            offset.push_back(pos);
//                        }
//                    }
//                }
//            }
//        }
        offset.resize(0);
        for(int i = 0; i < params->get_num_particles(); i++){
            dVector pos(0);
            for(int ndim = 0; ndim < params->get_ndim(); ndim++){
                pos.push_back(util->randnormed(params->get_box_size()));
            }
            offset.push_back(pos);
        }
    }
    
    beads = list_ptr(new PathList<dVector>(params->get_ndim(), params->get_box_size()));
    
    for(int slice = 0; slice < params->get_num_timeslices(); slice++){
        for(int ptcl = 0; ptcl < params->get_num_particles(); ptcl++){
            beads->push_back(offset[ptcl],ptcl);
        }
    }
    
    //beads->generate_separations();
    beads->generate_neighbors();
    
    
    if(params->is_charged()){
        int numcharges = params->get_num_chgs();
        for(int i = 0; i < params->get_num_particles();i++){
            charge_list.push_back(2*(i%numcharges)-1);
        }
    }
    
    //Make the beads a periodic chain
    beads->make_circular();
    worm = false;
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

void Path::put_worm_in_box(){
    
    std::vector<std::pair<int, int> > ht = beads->get_worm_indices();
    
    int start_row = ht[0].first;
    int start_slice = ht[0].second;
    int end_row = ht[1].first;
    int end_slice = ht[1].second;
    
    int slice = start_slice;
    int row = start_row;
    
    if(params->get_box_size() != -1){
        while(!(row == end_row && slice == end_slice)){
            beads->set_worm_bead_data(row, slice, util->location(beads->get_worm_bead_data(row, slice), params->get_box_size()),beads->get_worm_bead_data(row, slice));
            slice++;
            if(slice >= params->get_num_timeslices()){
                slice = slice%params->get_num_timeslices();
                row++;
            }
        }
    }
    
    
}

void Path::put_worm_in_box(iVector changed_rows, std::vector<std::pair<int, int> > start_end){
        
    if(params->get_box_size() != -1){
        for(iVector::iterator row = changed_rows.begin(); row != changed_rows.end(); row++){
            for(int slice = start_end[row-changed_rows.begin()].first; slice <= start_end[row-changed_rows.begin()].second; slice++){
                beads->set_worm_bead_data(*row, slice, util->location(beads->get_worm_bead_data(*row, slice), params->get_box_size()),beads->get_worm_bead_data(*row, slice));
            }
        }
    }
    

}