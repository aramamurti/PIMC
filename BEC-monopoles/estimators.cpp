//
//  estimators.cpp
//  BEC-monopoles
//
//  Created by Adith Ramamurti on 11/12/15.
//  Copyright © 2015 Adith Ramamurti. All rights reserved.
//

#include "estimators.hpp"

float Energy_Estimator::total_energy(){
    float ke =  kinetic_energy();
    float pe = potential_energy();
    
    return ke +pe;
}

float Energy_Estimator::potential_energy(){
    float pe = 0;
    int num_particles = path->get_beads()->get_num_particles();
    for(std::vector<boost::shared_ptr<Potential_Functions> >::iterator it = pot_funcs.begin(); it != pot_funcs.end(); it++){
        switch(potentials[it-pot_funcs.begin()]){
            case 0:
                for(int slice = 0; slice < num_timeslices; slice ++)
                    for(int ptcl = 0; ptcl < num_particles; ptcl ++)
                        pe += (*it)->potential_value(path->get_beads()->get_bead_data(ptcl, slice));
                break;
            case 1:
            case 2:
                for(int slice = 0; slice < num_timeslices; slice ++)
                    for(int ptcl = 0; ptcl < num_particles; ptcl ++)
                        for(int i = ptcl+1; i < num_particles; i++){
                            fVector distvec =utility->dist(path->get_beads()->get_pair_same_slice(ptcl, i, slice), path->get_parameters()->get_box_size());
                            float dist = sqrt(inner_product(distvec.begin(), distvec.end(),distvec.begin(), 0.0));
                            pe += (*it)->potential_value(dist);
                        }
                break;
            case 3:
                break;
        }
    }
    
    return pe/num_timeslices;
    
}

float Energy_Estimator::kinetic_energy(){
    int num_particles = path->get_beads()->get_num_particles();
    float ke = 0.0;
    float tot = 0.0;
    for(int slice = 0; slice < num_timeslices; slice++){
        for(int ptcl = 0; ptcl < num_particles; ptcl++){
            ffVector pair = path->get_beads()->get_pair_same_path(ptcl, slice, 1);
            fVector distVec = utility->dist(pair, path->get_parameters()->get_box_size());
            float dist = inner_product(distVec.begin(), distVec.end(), distVec.begin(), 0.0);
            tot -= norm*dist;
        }
    }
    
    ke = 0.5*path->get_parameters()->get_ndim()*num_particles/path->get_parameters()->get_tau() +tot/num_timeslices;
    
    return ke;
}