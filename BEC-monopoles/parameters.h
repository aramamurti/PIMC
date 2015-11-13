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
    float tau;
    float T;
    float kb;
    float lambda;
    int ndim;
    int timeslices;
    int particles;
    bool boson;
    int end_step;
    int skip;
    int equilibration;
    float box_size;
    bool per_bound_cond;
    std::vector<bool> potentials;
    int charges;
    bool charged;

    
public:
    Parameters(){
        ndim = 1;
        kb = 1.0;
        T = 1.0;
        lambda = 0.5;//pow(hbar,2)/(2*m);
        
        charges = 1;
        charged = true;
        
        boson = true;
        
        particles = 4;
        timeslices = 10;
        end_step = 10000;
        skip = 10;
        equilibration = 1000;
        
        potentials.resize(4);
        potentials[0] = true;
        
        tau = 1/(T*timeslices);
        
        per_bound_cond = false;
        if(per_bound_cond)
            box_size = pow(particles/10.,1/3.)*7.7099;
        else
            box_size = -1;

    }
    ~Parameters(){}
    
    void set_T(float newT){
        T = newT;
        tau = 1/(T*timeslices);
    }
    void set_timeslices(float newTS){
        timeslices = newTS;
        tau = 1/(T*timeslices);
    }
    
    int get_ndim(){return ndim;}
    float get_T(){return T;}
    float get_kb(){return kb;}
    bool is_boson(){return boson;}
    bool is_charged(){return charged;}
    int get_num_chgs(){return charges;}
    float get_tau(){return tau;}
    float get_lambda(){return lambda;}
    int get_end_step(){return end_step;}
    int get_skip(){return skip;}
    int get_equilibration(){return equilibration;}
    int get_num_timeslices(){return timeslices;}
    int get_num_particles(){return particles;}
    float get_box_size(){return box_size;}
    std::vector<bool> get_potentials(){return potentials;}
};

#endif /* defined(__BEC_monopoles__parameters__) */
