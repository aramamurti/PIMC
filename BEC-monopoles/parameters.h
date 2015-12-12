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

/*************************************************************************************
 
                                PARAMETERS CLASS
 
 Class Description: This class takes parameters from a parameters map and sets the
    appropriate values for the simulation (moves, potentials, box size, etc.). The
    class also handles the getting of parameters by other objects in the simulation.
 
 *************************************************************************************/

class Parameters{
    
private:
    
    //Map of parameters that are read in from file. Key is name of parameter, and value is parameter value
    typedef boost::unordered_map<std::string, std::string> sectionmap;
    
    //Simulation parameters
    double tau, T, kb, lambda, box_size, mu, C0;
    int ndim, timeslices, particles, skip, equilibration, end_step, charges, mbar;
    bool boson, per_bound_cond, charged;
    
    std::string particle_type;
    
    std::vector<bool> potentials, move_list;

    
public:
    
    /*****************************************************************************************

                                        CONSTRUCTOR
     
     Method Description: This method takes a parameters map, finds the appropriate parameter 
        values and sets the simulation parameters accordingly. If the parameter isn't found in
        the map, a default value is set.
     
     *****************************************************************************************/
    Parameters(boost::unordered_map<std::string, boost::shared_ptr<sectionmap> > parameter_map){
        sectionmap consts = *(*parameter_map.find("Constants")).second;
        sectionmap moves = *(*parameter_map.find("Moves")).second;
        sectionmap simpars = *(*parameter_map.find("SimulationParameters")).second;
        
        sectionmap::iterator param_it;
        
        //Number of dimensions
        param_it = consts.find("ndim");
        if(param_it != consts.end())
            ndim = std::stoi((*param_it).second);
        else
            ndim = 1;
        
        //Boltzmann constant
        param_it = consts.find("kb");
        if(param_it != consts.end())
            kb = std::stod((*param_it).second);
        else
            kb = 1;
        
        //Worm Constant
        param_it = consts.find("C0");
        if(param_it != consts.end())
            C0 = std::stod((*param_it).second);
        else
            C0 = 0;
        
        //Number of particles
        param_it = simpars.find("particles");
        if(param_it != simpars.end())
            particles = std::stoi((*param_it).second);
        else
            particles = 1;
        
        
        //Equilibration steps
        param_it = simpars.find("equilibration");
        if(param_it != simpars.end())
            equilibration = std::stoi((*param_it).second);
        else
            equilibration = 1000;
        
        
        //Simulation end step
        param_it = simpars.find("end_step");
        if(param_it != simpars.end())
            end_step = std::stoi((*param_it).second);
        else
            end_step = 10000;
        
        //Temperature
        param_it = simpars.find("temperature");
        if(param_it != simpars.end())
            T = std::stod((*param_it).second);
        else
            T = 1.0;
        
        //Time slices
        param_it = simpars.find("time_slices");
        if(param_it != simpars.end())
            set_timeslices(std::stoi((*param_it).second));
        else
            set_timeslices(100);
        
        //Chemical Potential
        param_it = simpars.find("mu");
        if(param_it != simpars.end())
            mu = std::stod((*param_it).second);
        else
            mu = 0;
        
        
        //Periodic boundary conditions
        param_it = simpars.find("periodic_boundary_conditions");
        per_bound_cond = false;
        if(param_it != simpars.end())
            per_bound_cond = ((*param_it).second == "true");
        
        //Number of different charges
        charges = 0;
        
        //Default box size
        box_size = -1;
        
        //Possible potentials
        potentials.resize(4);

        //Particle properties set by type
        param_it = simpars.find("particle_type");
        if(param_it != simpars.end()){
            if((*param_it).second == "boson_harmonic"){
                particle_type = (*param_it).second;
                lambda = 0.5;
                charged = false;
                boson = true;
                if(per_bound_cond)
                    box_size = 10;
                potentials[0] = true;

            }
            else if((*param_it).second == "boltzmannon_harmonic"){
                particle_type = (*param_it).second;
                lambda = 0.5;
                charged = false;
                boson = false;
                if(per_bound_cond)
                    box_size = 100;
                potentials[0] = true;
            }
            else if((*param_it).second == "he4"){
                particle_type = (*param_it).second;
                lambda = 6.0596;
                charged = false;
                boson = true;
                if(per_bound_cond)
                    box_size = pow(particles/10.,1./ndim)*7.7099;
                potentials[2] = true;
            }
            else if((*param_it).second == "boson_coulomb"){
                particle_type = (*param_it).second;
                lambda = 0.5;
                charged = true;
                charges = 1;
                boson = true;
                if(per_bound_cond)
                    box_size = (particles)^(1/ndim);
                potentials[3] = true;
            }
            else if((*param_it).second == "monopole_liquid"){
                particle_type = (*param_it).second;
                lambda = 0.5;
                charged = true;
                charges = 2;
                boson = true;
                if(per_bound_cond)
                    box_size = (particles)^(1/ndim);
                potentials[3] = true;
            }
                
        }
        
        //Possible moves allowed
        
        param_it = moves.find("Center_of_Mass");
        if(param_it != moves.end())
            move_list.push_back((*param_it).second == "true");
        
        param_it = moves.find("Bisection");
        if(param_it != moves.end())
            move_list.push_back((*param_it).second == "true");
        
        param_it = moves.find("Worm");
        if(param_it != moves.end()){
            if(!boson){
                move_list.push_back(false);
                move_list.push_back(false);
            }
            else if((*param_it).second == "true"){
                move_list.push_back(false);
                move_list.push_back(true);
            }
            else{
                move_list.push_back(true);
                move_list.push_back(false);
            }
        }
        
        
    }
    
    /*************************************************************************
     
                                DESTRUCTOR
     
     *************************************************************************/
    
    ~Parameters(){}
    
    /*************************************************************************
     
     Method Description: This method adjusts the number of timeslices for each
        particle in the system, and adjusts tau accordingly.
     
     *************************************************************************/
    
    void set_timeslices(int new_ts){
        timeslices = new_ts;
        tau = 1.0/(T*timeslices);        
        int mbar_try = timeslices/2;
        if(mbar_try%2)
            mbar_try++;
        set_mbar(mbar_try);
            
    }
    
    void shift_timeslices(int shift){
        set_timeslices(timeslices + shift);
    }
    
    void set_mbar(int new_mbar){
        mbar = new_mbar;
    }
    
    void shift_mu(double shift){
        mu += shift;
    }
    
    void shift_C0(double shift){
        C0 += shift*C0;
    }
    
    /*************************************************************************
     
                                GETTER METHODS
     
     *************************************************************************/
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
    double get_mu(){return mu;}
    double get_C0(){return C0;}
    int get_mbar(){return mbar;}
    std::string get_particle_type(){return particle_type;}
    std::vector<bool> get_potentials(){return potentials;}
    std::vector<bool> get_move_list(){return move_list;}
};

#endif /* defined(__BEC_monopoles__parameters__) */
