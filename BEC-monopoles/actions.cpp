//
//  actions.cpp
//  BEC-monopoles
//
//  Created by Adith Ramamurti on 11/12/15.
//  Copyright © 2015 Adith Ramamurti. All rights reserved.
//

#include "actions.hpp"

double Potential_Action::potential_helper(int slice, int ptcl){
    double pe = 0;
    int num_particles = path->get_beads()->get_num_particles();
    
    for(std::vector<boost::shared_ptr<Potential_Functions> >::iterator it = pot_funcs.begin(); it != pot_funcs.end(); it++){
        switch(potentials[it-pot_funcs.begin()]){
            case 0:
                pe += (*it)->potential_value(path->get_beads()->get_bead_data(ptcl, slice));
                break;
            case 1:
            case 2:
                for(int i = ptcl+1; i < num_particles; i++){
                    const dVector& distvec = path->get_beads()->get_path_separation(ptcl, i, slice);
                    double dist = sqrt(inner_product(distvec.begin(), distvec.end(),distvec.begin(), 0.0));
                    pe += (*it)->potential_value(dist);
                }
                break;
            case 3:
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
    return pe;
}


double Potential_Action::get_action(int slice, int dist, int start_omit, int end_omit){
    double pot = 0;
    int num_particles = path->get_beads()->get_num_particles();
    for(int ptcl = 0; ptcl < num_particles; ptcl++){
        pot += potential_helper(slice, ptcl);
    }
    
    return path->get_parameters()->get_tau()*pot;
}

double Potential_Action::get_action_single_particle(int ptcl, int slice){
    double pot = 0;
    int num_particles = path->get_beads()->get_num_particles();
    for(std::vector<boost::shared_ptr<Potential_Functions> >::iterator it = pot_funcs.begin(); it != pot_funcs.end(); it++){
        switch(potentials[it-pot_funcs.begin()]){
            case 0:
                pot += (*it)->potential_value(path->get_beads()->get_bead_data(ptcl, slice));
                break;
            case 1:
            case 2:
                for(int i = 0; i < num_particles; i++)
                    if(i != ptcl){
                        const dVector& distvec = path->get_beads()->get_path_separation(ptcl, i, slice);
                        double dist = sqrt(inner_product(distvec.begin(), distvec.end(),distvec.begin(), 0.0));
                        pot += (*it)->potential_value(dist);
                    }
                break;
            case 3:
                for(int i = 0; i < num_particles; i++){
                    const dVector& distvec = path->get_beads()->get_path_separation(ptcl, i, slice);
                    int chgi = path->get_beads()->get_charge(ptcl);
                    int chgj = path->get_beads()->get_charge(i);
                    pot += (*it)->potential_value(distvec, chgi, chgj, path->get_parameters()->get_box_size());
                }
        }
    }
    return path->get_parameters()->get_tau()*pot;
}

double Potential_Action::get_action_multiple_particles(iVector ptcls, int slice){
    double pot = 0;
    int num_particles = path->get_beads()->get_num_particles();
    for(std::vector<boost::shared_ptr<Potential_Functions> >::iterator it = pot_funcs.begin(); it != pot_funcs.end(); it++){
        switch(potentials[it-pot_funcs.begin()]){
            case 0:
                for(iVector::iterator ptcl = ptcls.begin(); ptcl != ptcls.end(); ptcl++){
                    pot += (*it)->potential_value(path->get_beads()->get_bead_data(*ptcl, slice));
                }
                break;
            case 1:
            case 2:
                for(iVector::iterator ptcl = ptcls.begin(); ptcl != ptcls.end(); ptcl++){
                    for(int i = 0; i < num_particles; i++)
                        if(i != *ptcl){
                            const dVector& distvec = path->get_beads()->get_path_separation(*ptcl, i, slice);
                            double dist = sqrt(inner_product(distvec.begin(), distvec.end(),distvec.begin(), 0.0));
                            pot += (*it)->potential_value(dist);
                        }
                }
                for(iVector::iterator ptcl = ptcls.begin(); ptcl != ptcls.end(); ptcl++){
                    for(iVector::iterator ptcl2 = ptcl+1; ptcl2 != ptcls.end(); ptcl2++){
                        const dVector& distvec = path->get_beads()->get_path_separation(*ptcl, *ptcl2, slice);
                        double dist = sqrt(inner_product(distvec.begin(), distvec.end(),distvec.begin(), 0.0));
                        pot -= (*it)->potential_value(dist);
                    }
                }
                break;
            case 3:
                for(iVector::iterator ptcl = ptcls.begin(); ptcl != ptcls.end(); ptcl++){
                    for(int i = 0; i < num_particles; i++){
                        const dVector& distvec = path->get_beads()->get_path_separation(*ptcl, i, slice);
                        int chgi = path->get_beads()->get_charge(*ptcl);
                        int chgj = path->get_beads()->get_charge(i);
                        pot += (*it)->potential_value(distvec, chgi, chgj, path->get_parameters()->get_box_size());
                    }
                }
                for(iVector::iterator ptcl = ptcls.begin(); ptcl != ptcls.end(); ptcl++){
                    for(iVector::iterator ptcl2 = ptcl+1; ptcl2 != ptcls.end(); ptcl2++){
                        const dVector& distvec = path->get_beads()->get_path_separation(*ptcl, *ptcl2, slice);
                        int chgi = path->get_beads()->get_charge(*ptcl);
                        int chgj = path->get_beads()->get_charge(*ptcl2);
                        pot -= (*it)->potential_value(distvec, chgi, chgj, path->get_parameters()->get_box_size());
                    }
                }
                break;
        }
    }
    return path->get_parameters()->get_tau()*pot;
}