//
//  parameters.h
//  BEC-monopoles
//
//  Created by Adith Ramamurti on 5/29/15.
//  Copyright (c) 2015 Adith Ramamurti. All rights reserved.
//

#ifndef __BEC_monopoles__parameters__
#define __BEC_monopoles__parameters__

#include "uni_header.h"

class Parameters{
    
private:
    
    typedef boost::unordered_map<std::string, std::string> sectionmap;
    
    double tau, T, kb, lambda, box_size;
    int ndim, timeslices, particles, skip, equilibration, end_step, charges;
    bool boson, per_bound_cond, charged;
    std::vector<bool> potentials, move_list;
    std::string particle_type;

    
public:
    
    Parameters(boost::unordered_map<std::string, boost::shared_ptr<sectionmap> > parameter_map){
        sectionmap consts = *(*parameter_map.find("Constants")).second;
        sectionmap moves = *(*parameter_map.find("Moves")).second;
        sectionmap simpars = *(*parameter_map.find("SimulationParameters")).second;
        
        sectionmap::iterator param_it;
        
        param_it = consts.find("ndim");
        if(param_it != consts.end())
            ndim = std::stoi((*param_it).second);
        else
            ndim = 1;
        param_it = consts.find("kb");
        if(param_it != consts.end())
            kb = std::stof((*param_it).second);
        else
            kb = 1;
        
        param_it = moves.find("Center_of_Mass");
        if(param_it != moves.end())
            move_list.push_back((*param_it).second == "true");
        
        param_it = moves.find("Bisection");
        if(param_it != moves.end())
            move_list.push_back((*param_it).second == "true");

        param_it = moves.find("Perm_Bisection");
        if(param_it != moves.end())
            move_list.push_back((*param_it).second == "true");

        
        param_it = simpars.find("particles");
        if(param_it != simpars.end())
            particles = std::stoi((*param_it).second);
        else
            particles = 1;
        
        param_it = simpars.find("equilibration");
        if(param_it != simpars.end())
            equilibration = std::stoi((*param_it).second);
        else
            equilibration = 1000;
        
        param_it = simpars.find("end_step");
        if(param_it != simpars.end())
            end_step = std::stoi((*param_it).second);
        else
            end_step = 10000;
        
        param_it = simpars.find("temperature");
        if(param_it != simpars.end())
            T = std::stof((*param_it).second);
        else
            T = 1.0;
        
        param_it = simpars.find("periodic_boundary_conditions");
        per_bound_cond = false;
        if(param_it != simpars.end())
            per_bound_cond = ((*param_it).second == "true");
        
        charges = 0;
        box_size = -1;
        potentials.resize(4);

        param_it = simpars.find("particle_type");
        if(param_it != simpars.end()){
            if((*param_it).second == "boson_harmonic"){
                lambda = 0.5;
                charged = false;
                boson = true;
                if(per_bound_cond)
                    box_size = 10;
                potentials[0] = true;

            }
            else if((*param_it).second == "boltzmannon_harmonic"){
                lambda = 0.5;
                charged = false;
                boson = false;
                if(per_bound_cond)
                    box_size = 10;
                potentials[0] = true;
            }
            else if((*param_it).second == "he4"){
                lambda = 6.0596;
                charged = false;
                boson = true;
                if(per_bound_cond)
                    box_size = pow(particles/10.,1./ndim)*7.7099;
                potentials[2] = true;
            }
            else if((*param_it).second == "boson_coulomb"){
                lambda = 0.5;
                charged = true;
                charges = 1;
                boson = true;
                if(per_bound_cond)
                    box_size = (particles)^(1/ndim);
                potentials[3] = true;
            }
            else if((*param_it).second == "monopole_liquid"){
                lambda = 0.5;
                charged = true;
                charges = 2;
                boson = true;
                if(per_bound_cond)
                    box_size = (particles)^(1/ndim);
                potentials[3] = true;
            }
                
        }
        
        set_timeslices(40);
        skip = 10;
        
    }
    
    ~Parameters(){}
    
    void set_timeslices(double newTS){
        timeslices = newTS;
        tau = 1/(T*timeslices);
    }
    
    int get_ndim(){return ndim;}
    double get_T(){return T;}
    double get_kb(){return kb;}
    bool is_boson(){return boson;}
    bool is_charged(){return charged;}
    int get_num_chgs(){return charges;}
    double get_tau(){return tau;}
    double get_lambda(){return lambda;}
    int get_end_step(){return end_step;}
    int get_skip(){return skip;}
    int get_equilibration(){return equilibration;}
    int get_num_timeslices(){return timeslices;}
    int get_num_particles(){return particles;}
    double get_box_size(){return box_size;}
    std::vector<bool> get_potentials(){return potentials;}
    std::vector<bool> get_move_list(){return move_list;}
};

#endif /* defined(__BEC_monopoles__parameters__) */
