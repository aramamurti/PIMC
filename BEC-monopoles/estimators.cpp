//
//  estimators.cpp
//  BEC-monopoles
//
//  Created by Adith Ramamurti on 11/12/15.
//  Copyright © 2015 Adith Ramamurti. All rights reserved.
//

#include "estimators.hpp"

fVector Energy_Estimator::estimate(){
    float ke =  kinetic_energy();
    float pe = potential_energy();
    
    fVector energy;
    energy.push_back(ke+pe);
    energy.push_back(ke);
    energy.push_back(pe);
    
    return energy;
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

fVector Winding_Estimator::estimate(){
    int num_particles = path->get_beads()->get_num_particles();
    fVector dvectot(ndim,0.0);
    for(int ptcl = 0; ptcl < num_particles; ptcl++){
        for(int slice = 0; slice < num_timeslices; slice++){
            ffVector pair = path->get_beads()->get_pair_same_path(ptcl, slice, 1);
            fVector distVec = utility->dist(pair, path->get_parameters()->get_box_size());
            dvectot = utility->vecadd(dvectot, distVec);
        }
    }
    iVector wnum(ndim,0);
    for(fVector::iterator it = dvectot.begin(); it != dvectot.end(); it++){
        wnum[it-dvectot.begin()] = (int)round(*it/path->get_parameters()->get_box_size());
    }
    return fVector(wnum.begin(),wnum.end());
}

fVector Permutation_Estimator::estimate(){
    iVector cycles = path->get_beads()->get_cycles();
    return fVector(cycles.begin(), cycles.end());
}