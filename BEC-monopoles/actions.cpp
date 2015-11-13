//
//  actions.cpp
//  BEC-monopoles
//
//  Created by Adith Ramamurti on 11/12/15.
//  Copyright © 2015 Adith Ramamurti. All rights reserved.
//

#include "actions.hpp"

float Potential_Action::potential_helper(int slice, int ptcl){
    float pe = 0;
    int num_particles = path->get_beads()->get_num_particles();
    
    for(std::vector<boost::shared_ptr<Potential_Functions> >::iterator it = pot_funcs.begin(); it != pot_funcs.end(); it++){
        switch(potentials[it-pot_funcs.begin()]){
            case 0:
                pe += 2*(*it)->potential_value(path->get_beads()->get_bead_data(ptcl, slice));
                break;
            case 1:
            case 2:
                for(int i = 0; i < num_particles; i++)
                    if(i!=ptcl){
                        fVector distvec =utility->dist(path->get_beads()->get_pair_same_slice(ptcl, i, slice), path->get_parameters()->get_box_size());
                        float dist = sqrt(inner_product(distvec.begin(), distvec.end(),distvec.begin(), 0.0));
                        pe += (*it)->potential_value(dist);
                    }
                break;
            case 3:
                break;
        }
    }
    
    pe = 0.5*pe;
    return pe;
}

float Potential_Action::get_action(int slice, int dist){
    float pot = 0;
    int num_particles = path->get_beads()->get_num_particles();
    for(int ptcl = 0; ptcl < num_particles; ptcl++){
        pot += potential_helper(slice, ptcl);
    }
    return path->get_parameters()->get_tau()*pot;
}

float Kinetic_Action::get_action(int slice, int dist){
    float kin = 0;
    int num_particles = path->get_beads()->get_num_particles();

    for(int ptcl = 0; ptcl < num_particles; ptcl++){
        ffVector pair = path->get_beads()->get_pair_same_path(ptcl, slice, dist);
        fVector distVec = utility->dist(pair, path->get_parameters()->get_box_size());
        float ipdist =  inner_product(distVec.begin(),distVec.end(),distVec.begin(),0.0);
        kin += norm/dist*ipdist;
    }
    return kin;
}