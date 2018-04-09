//
//  parameters.hpp
//
//
//  Created by Adith Ramamurti on 9/12/16.
//
//

#ifndef parameters_h
#define parameters_h

#include <cmath>
#include <boost/any.hpp>
#include <boost/unordered_map.hpp>


/*---------------------------------------------------------------------------------------------------*

This file contains the Parameters class, which handles storing, as well as reading in, the parameters
of the simulation.

*---------------------------------------------------------------------------------------------------*/

typedef boost::unordered_map<std::string, std::string> sectionmap;

class Parameters{
public:
    int total_slices = -1;
    int my_start = -1; //processor start slice
    int my_end = -1; //processor end slice
    int slices_per_process = -1;
    int num_workers = 0; //number of processors
    int multistep_dist = 2;
    
    int equilibration = 0;
    int end = 0;
    
    int dimensions = 1;
    int particles = 1;
    int init_particles = 1;
    double real_particles = 1;
    double lambda = 1;
    double tau = .1;
    double temperature = 1;
    double kb = 1;
    
    bool per_bound_cond;
    double box_size = -1;
    double volume = -1;
    double volume2 = -1;
    double grid_size = 0;
    
    std::string particle_type;
    bool boson = false;
    
    bool charged = false;
    int charges = 0;
    double coupling = 0;
    
    //Worm parameters
    double C0 = 1;
    double C = 1;
    int Mbar = 2;
    double mu = 0;
    bool gce = false;
    bool worm_on = false;
    std::pair<int, int> worm_head; //row, col of worm head
    std::pair<int, int> worm_tail; //row, col of worm tail
    int worm_length = 0;
    
    int potential = 0;
    
    //Ewald parameters
    double p = 23.0258;//-log(ewaldError);
    double alpha = -1;
    double coulcut = -1;
    double coulcut2 = -1;
    double kcut = -1;
    double kcut2 = -1;
    int nmax = -1;
    int kmax = -1;
    double kfac = -1;
    
    Parameters(){
        total_slices = 0;
        particles = 0;
    }
    
    void set_parameters(boost::unordered_map<std::string, sectionmap> parameter_map){
        sectionmap simpars = (*parameter_map.find("SimulationParameters")).second;
        
        sectionmap::iterator param_it;
        
        //Number of dimensions
        param_it = simpars.find("dimensions");
        if(param_it != simpars.end())
            dimensions = std::stoi((*param_it).second);
        else
            dimensions = 1;
        
        //Boltzmann constant
        param_it = simpars.find("kb");
        if(param_it != simpars.end())
            kb = std::stod((*param_it).second);
        else
            kb = 1;
        
        
        //Number of particles
        param_it = simpars.find("particles");
        if(param_it != simpars.end())
            particles = std::stoi((*param_it).second);
        else
            particles = 1;
        
        init_particles = particles;
        real_particles = particles;
        
        //Equilibration steps
        param_it = simpars.find("equilibration");
        if(param_it != simpars.end())
            equilibration = std::stoi((*param_it).second);
        else
            equilibration = 1000;
        
        
        //Simulation end step
        param_it = simpars.find("end_step");
        if(param_it != simpars.end())
            end = std::stoi((*param_it).second);
        else
            end = 10000;
        
        //Temperature
        param_it = simpars.find("temperature");
        if(param_it != simpars.end())
            temperature = std::stod((*param_it).second);
        else
            temperature = 1.0;
        
        //Time slices
        param_it = simpars.find("time_slices");
        if(param_it != simpars.end())
            set_total_slices(std::stoi((*param_it).second));
        else
            set_total_slices(100);
        
        //Periodic boundary conditions
        param_it = simpars.find("periodic_boundary_conditions");
        per_bound_cond = false;
        if(param_it != simpars.end())
            per_bound_cond = ((*param_it).second == "true");
        
        //Temperature
        param_it = simpars.find("NN_grid_size");
        if(param_it != simpars.end())
            grid_size = std::stod((*param_it).second);
        else
            grid_size = 1.0;

        //Coupling constant
        param_it = simpars.find("coupling");
        if(param_it != simpars.end())
            coupling = std::stod((*param_it).second);
        else
            coupling = 0;
        
        //Chemical potential
        param_it = simpars.find("mu");
        if(param_it != simpars.end())
            mu = std::stod((*param_it).second);
        else
            mu = 0;
        
        //Worm constant
        param_it = simpars.find("C0");
        if(param_it != simpars.end())
            set_C0(std::stod((*param_it).second));
        else
            set_C0(0);
        
        //Worm
        param_it = simpars.find("grand_canonical_ensemble");
        gce = false;
        if(param_it != simpars.end())
            gce = ((*param_it).second == "true");
        
        //Number of different charges
        charges = 0;
        
        //Default box size
        box_size = -1;
        
        
        //Particle properties set by type
        param_it = simpars.find("particle_type");
        if(param_it != simpars.end()){
            if((*param_it).second == "boson_harmonic"){
                particle_type = (*param_it).second;
                lambda = 0.5;
                charged = false;
                boson = true;
                if(per_bound_cond)
                    box_size = 30;
                potential = 0;
            }
            else if((*param_it).second == "boltzmannon_harmonic"){
                particle_type = (*param_it).second;
                lambda = 0.5;
                charged = false;
                boson = false;
                if(per_bound_cond)
                    box_size = 30;
                potential = 0;
            }
            else if((*param_it).second == "he4"){
                particle_type = (*param_it).second;
                lambda = 6.0596;
                charged = false;
                boson = true;
                if(per_bound_cond){
                    box_size = pow(particles/10.,1./dimensions)*7.7099;
                    volume = pow(box_size, 3);
                    volume2 = pow(volume, 2);
                    kfac = 2*M_PI/box_size;
                }
                potential = 1;
            }
            else if((*param_it).second == "boson_coulomb"){
                particle_type = (*param_it).second;
                lambda = 0.5;
                charged = true;
                charges = 1;
                boson = true;
                if(per_bound_cond){
                    box_size = pow(particles, 1./dimensions);
                    volume = pow(box_size, 3);
                    volume2 = pow(volume, 2);
                    kfac = 2*M_PI/box_size;
                }
                potential = 2;
            }
            else if((*param_it).second == "monopole_liquid"){
                particle_type = (*param_it).second;
                lambda = 0.5;
                charged = true;
                charges = 2;
                boson = true;
                if(per_bound_cond){
                    box_size = pow(particles/2.,1./dimensions);
                    volume = pow(box_size, 3);
                    volume2 = pow(volume, 2);
                    kfac = 2*M_PI/box_size;
                }
                potential = 2;
            }	
            else if((*param_it).second == "scaled_monopole_liquid"){
                particle_type = (*param_it).second;
                lambda = 0.5;
                charged = true;
                charges = 2;
                boson = true;
                if(per_bound_cond){
                    box_size = std::pow((1./(3.191667129161336*0.557/std::pow(std::log(temperature/3.45*2.69),2.)*std::pow(temperature/3.45,3.))*particles/2.),1./3.);
                    volume = pow(box_size, 3);
                    volume2 = pow(volume, 2);
                    kfac = 2*M_PI/box_size;
                }
                potential = 2;
            }
        }
    }
    
    void set_total_slices(int new_ts){
        total_slices = new_ts;
        tau = 1.0/(temperature*total_slices);
        int Mbar_try = total_slices/5;
        if(Mbar_try%2)
            ++Mbar_try;
        Mbar = Mbar_try;
    }
    
    void set_temperature(double T){
        temperature = T;
        tau = 1.0/(temperature*total_slices);
    }
    
    void set_C(double newC){
        C = newC;
        if(box_size != -1)
            C0 = C*pow(box_size,3)*total_slices*Mbar;
        else
            C0 = C*pow(2,3)*total_slices*Mbar;
    }
    
    void set_C0(double newC0){
        C0 = newC0;
        if(box_size != -1)
            C = C0/(pow(box_size,3)*total_slices*Mbar);
        else
            C = C0/(pow(2,3)*total_slices*Mbar);
    }
    
    void shift_mu(double shift){
        mu+=shift;
    }
    
    void shift_C0(double shift){
        set_C0(C0+C0*shift);
    }
    
    void set_multistep(){
        multistep_dist = pow(2, floor(log2(total_slices/4.)));
        if(multistep_dist < 2) multistep_dist = 2;
    }
    
};

#endif /* parameters_h */
