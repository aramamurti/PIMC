//
//  estimators.cpp
//  BEC-monopoles
//
//  Created by Adith Ramamurti on 11/12/15.
//  Copyright © 2015 Adith Ramamurti. All rights reserved.
//

#include "estimators.hpp"

dVector Energy_Estimator::estimate(){
    double ke =  kinetic_energy();
    double pe = potential_energy();
    
    dVector energy;
    energy.push_back(ke+pe);
    energy.push_back(ke);
    energy.push_back(pe);
    
    return energy;
}

double Energy_Estimator::potential_energy(){
    double pe = 0;
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
                            const dVector& distvec = path->get_beads()->get_path_separation(ptcl, i, slice);
                            double dist = sqrt(inner_product(distvec.begin(), distvec.end(),distvec.begin(), 0.0));
                            pe += (*it)->potential_value(dist);
                        }
                break;
            case 3:
                for(int slice = 0; slice < num_timeslices; slice ++)
                    for(int ptcl = 0; ptcl < num_particles; ptcl ++)
                        for(int i = 0; i < num_particles; i++){
                            const dVector& distvec = path->get_beads()->get_path_separation(ptcl, i, slice);
                            int chgi = path->get_beads()->get_charge(ptcl);
                            int chgj = path->get_beads()->get_charge(i);
                            pe += (*it)->potential_value(distvec, chgi, chgj, path->get_parameters()->get_box_size());
                        }
                pe = pe/2;
                break;
        }
    }
    
    return pe/num_timeslices;
    
}

double Energy_Estimator::kinetic_energy(){
    int num_particles = path->get_beads()->get_num_particles();
    double ke = 0.0;
    double tot = 0.0;
    for(int slice = 0; slice < num_timeslices; slice++){
        for(int ptcl = 0; ptcl < num_particles; ptcl++){
            dVector bead1;
            dVector bead2;
            path->get_beads()->get_pair_same_path(ptcl, slice, 1, bead1, bead2);
            dVector distVec;
            utility->dist(bead1, bead2, distVec, path->get_parameters()->get_box_size());
            double dist = inner_product(distVec.begin(), distVec.end(), distVec.begin(), 0.0);
            tot -= norm*dist;
        }
    }
    
    ke = 0.5*path->get_parameters()->get_ndim()*num_particles/path->get_parameters()->get_tau() +tot/num_timeslices;
    
    return ke;
}

dVector Winding_Estimator::estimate(){
    int num_particles = path->get_beads()->get_num_particles();
    dVector dvectot(ndim,0.0);
    for(int ptcl = 0; ptcl < num_particles; ptcl++){
        for(int slice = 0; slice < num_timeslices; slice++){
            dVector bead1;
            dVector bead2;
            path->get_beads()->get_pair_same_path(ptcl, slice, 1, bead1, bead2);
            dVector distVec;
            utility->dist(bead1, bead2, distVec, path->get_parameters()->get_box_size());
            dvectot = utility->vecadd(dvectot, distVec);
        }
    }
    iVector wnum(ndim,0);
    for(dVector::iterator it = dvectot.begin(); it != dvectot.end(); it++){
        wnum[it-dvectot.begin()] = (int)round(*it/path->get_parameters()->get_box_size());
    }
    return dVector(wnum.begin(),wnum.end());
}

dVector Permutation_Estimator::estimate(){
    iVector cycles = path->get_beads()->get_cycles();
    return dVector(cycles.begin(), cycles.end());
}