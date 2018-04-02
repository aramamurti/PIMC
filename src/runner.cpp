//
//  runner.cpp
//  PIMC-WORM
//
//  Created by Adith Ramamurti on 10/7/16.
//  Copyright © 2016 Adith Ramamurti. All rights reserved.
//
// This file contains the methods to run the simulation

#include "runner.hpp"

//Initializer
Runner::Runner(int &id, Parameters &params, MPI_Comm &local){
    local_comm = local;
    
    //based on parameters, choose the appropriate moves to include in the simulation
    
    //general moves
    moves.push_back(new Center_of_Mass(id, params, local));
    moves.push_back(new Bisection(id, params, local));
    
    if(params.charges == 2) // two-charge moves
        moves.push_back(new Pair_Center_of_Mass(id, params, local));
    
    if(!params.gce && params.boson) // non-worm boson moves
        moves.push_back(new Permutation_Bisection(id,params,local));
    
    if(params.gce && params.boson){ // worm boson  moves
        moves.push_back(new Open(id, params, local));
        moves.push_back(new Close(id, params, local));
        moves.push_back(new Insert(id, params, local));
        moves.push_back(new Remove(id, params, local));
        moves.push_back(new Advance_Tail(id, params, local));
        moves.push_back(new Advance_Head(id, params, local));
        moves.push_back(new Recede_Tail(id, params, local));
        moves.push_back(new Recede_Head(id, params, local));
        moves.push_back(new Swap_Tail(id, params, local));
        moves.push_back(new Swap_Head(id, params, local));
    }
    
    move_choices.resize(2); // move choices vector holds indices of moves in the move vector in two categories/types: worm moves and non-worm moves
    int index = 0;
    for(auto &move : moves){
        if(move.get_move_type().first)
            move_choices[0].push_back(index);
        if(move.get_move_type().second)
            move_choices[1].push_back(index);
        ++index;
    }
}

//This method equilibrates the simulation
void Runner::equilibrate(int &id, Parameters &params, Paths &paths, RNG &rng, Cos& cos){

    //set attempt and acceptance counters
    int com_att = 200;
    int com_acc = 0;
    int com_acc_pair = 0;
    int prev_num_attempts = 0;
    int prev_num_attempts_pair = 0;
    int off_diag_att = 3000;
    int off_diag_att_ctr = 0;
    int conf_att = 500;
    int conf_counter = 0;
    int diag_counter = 0;
    double cumulative_particles = 0;
    
    //checker flags used to check separation and neighbor tables, used only for debugging
    bool check1 = false;
    bool check2 = false;
    
    int initial_moves = 2;
    if(params.charges == 2)
        initial_moves = 3;
    
    //equilibration loop
    for(int step = 0; step < params.equilibration; ++step){
        if(id == 0  && (step+1)%100 == 0)
            std::cout << "Equilibration step " << step+1 << "/" <<params.equilibration << ":\t" << params.real_particles << " particles" << std::endl;
        //for the first third of equilibration, we only do center of mass and bisection moves
        if(double(step)/params.equilibration < 1/3.){
            for(int j = 0; j < params.particles; ++j){ // these two loops run n^2 moves per step (this can be changed)
                for(int n = 0; n < params.particles; ++n){
                    int success = 0;
                    do{
                        if(check1){
                            paths.check_nt();
                            paths.check_separations(params);
                            paths.check_separations_kinetic(params);
                            paths.check_charge();
                        }
                        int choice = rng.randint(initial_moves);
                        MPI_Bcast(&choice, 1, MPI_INT, 0, local_comm);
                        success = moves[choice].attempt(id, params, paths, rng, cos);
                    }while(success != 0);
                    if(moves[0].get_num_attempts()%com_att == 0 && moves[0].get_num_attempts() != prev_num_attempts){ //adjust center of mass delta based on acceptance ratio (aiming for 0.5)
                        com_acc = moves[0].get_num_accepts()-com_acc;
                        double com_acc_rat = double(com_acc)/com_att;
                        if (com_acc_rat < 0.2)
                            moves[0].shift_delta(-0.6);
                        else if (com_acc_rat < 0.3)
                            moves[0].shift_delta(-0.4);
                        else if (com_acc_rat < 0.4)
                            moves[0].shift_delta(-0.2);
                        else if(moves[0].get_delta() < params.box_size){
                            if (com_acc_rat > 0.6)
                                moves[0].shift_delta(0.2);
                            else if (com_acc_rat > 0.7)
                                moves[0].shift_delta(0.4);
                            else if (com_acc_rat > 0.8)
                                moves[0].shift_delta(0.6);
                        }
                        com_acc = moves[0].get_num_accepts();
                        prev_num_attempts = moves[0].get_num_attempts();
                    }
                    if(moves[0].get_delta() > params.box_size)
                        moves[0].set_delta(params.box_size);
                }
            }
        }
        //for the last 2/3 of equilibration, equilibrate all methods
        else {
            for(int j = 0; j < std::max(params.particles,params.init_particles); ++j){
                for(int n = 0; n < std::max(params.particles,params.init_particles); ++n){
                    int success = 0;
                    do{
                        if(check2){
                            paths.check_nt();
                            paths.check_separations(params);
                            paths.check_separations_kinetic(params);
                            paths.check_charge();
                        }
                        int choice = 0;
                        if(params.worm_on)
                            choice = move_choices[1][(rng.randint(move_choices[1].size()))];
                        else
                            choice = move_choices[0][(rng.randint(move_choices[0].size()))];
                        MPI_Bcast(&choice, 1, MPI_INT, 0, local_comm);
                        success = moves[choice].attempt(id, params, paths, rng, cos);
                    }while(success != 0);
                }
                if(params.charges == 2){ //adjust pair center of mass delta for 0.5 acceptance ratio
                    if(moves[2].get_num_attempts()%com_att == 0 && moves[2].get_num_attempts() != prev_num_attempts_pair){
                        com_acc_pair = moves[2].get_num_accepts()-com_acc_pair;
                        double com_acc_rat = double(com_acc_pair)/com_att;
                        if (com_acc_rat < 0.2)
                            moves[2].shift_delta(-0.6);
                        else if (com_acc_rat < 0.3)
                            moves[2].shift_delta(-0.4);
                        else if (com_acc_rat < 0.4)
                            moves[2].shift_delta(-0.2);
                        else if(moves[2].get_delta() < params.box_size){
                            if (com_acc_rat > 0.6)
                                moves[2].shift_delta(0.2);
                            else if (com_acc_rat > 0.7)
                                moves[2].shift_delta(0.4);
                            else if (com_acc_rat > 0.8)
                                moves[2].shift_delta(0.6);
                        }
                        com_acc_pair = moves[2].get_num_accepts();
                        prev_num_attempts_pair = moves[2].get_num_attempts();
                    }
                    if(moves[2].get_delta() > params.box_size)
                        moves[2].set_delta(params.box_size);
                }
                if(params.gce && params.boson){ //adjust worm constant and chemical potential to get appropriate number of diagonal configurations and particles in the system
                    params.real_particles = params.particles - (ceil(double(params.worm_length)/params.total_slices) - double(params.worm_length)/params.total_slices);
                    cumulative_particles += params.real_particles;
                    ++off_diag_att_ctr;
                    ++conf_counter;
                    if(!params.worm_on)
                        ++diag_counter;
                    if(off_diag_att_ctr == off_diag_att){
                        double avg_part = cumulative_particles/off_diag_att_ctr;
                        double mu_shift = 1 - avg_part/params.init_particles;
                        params.shift_mu(mu_shift);
                        cumulative_particles = 0;
                        off_diag_att_ctr = 0;
                    }
                    if(conf_counter == conf_att){
                        double diag_frac = double(diag_counter)/conf_counter;
                        if (params.C0 > 1.0E-5) {
                            if (diag_frac < 0.2)
                                params.shift_C0(-0.5);
                            else if (diag_frac >= 0.2 && diag_frac < 0.3)
                                params.shift_C0(-0.4);
                            else if (diag_frac >= 0.3 && diag_frac < 0.4)
                                params.shift_C0(-0.3);
                            else if (diag_frac >= 0.4 && diag_frac < 0.5)
                                params.shift_C0(-0.2);
                            else if (diag_frac >= 0.5 && diag_frac < 0.6)
                                params.shift_C0(-0.1);
                            else if (diag_frac >= 0.6 && diag_frac < 0.75)
                                params.shift_C0(-0.05);
                        }
                        if (params.C0 < 1.0E4) {
                            if (diag_frac <= 0.9 && diag_frac > 0.85)
                                params.shift_C0(0.05);
                            else if (diag_frac <= 0.95 && diag_frac > 0.9)
                                params.shift_C0(0.1);
                            else if (diag_frac > 0.95)
                                params.shift_C0(0.2);
                        }
                        conf_counter = 0;
                        diag_counter = 0;
                    }
                }
            }
        }
    }
}

//runs the simulation and writes output to files
void Runner::run(int &id, Parameters &params, IO &writer){

    //initialize path configurations and all helper methods/tables
    Paths paths(id, params, local_comm);
    RNG rng;
    rng.seed(id);
    Cos cos;
    cos.set_up();
    
    //initialize potential table for non-grand canonical simulations
    if(!params.gce || params.potential != 2)
        moves[0].initialize_potential_table(id, params, paths, cos);
    
    //equilibrate system
    if(id == 0)
        std::cout << "Equilibrating ... " << std::endl;
    equilibrate(id, params, paths, rng, cos);
    
    //write equilibration parameters
    if(id == 0){
        if(params.charges == 2)
            writer.write_equil_parameters(params, moves[0].get_delta(), moves[2].get_delta());
        else
            writer.write_equil_parameters(params, moves[0].get_delta());
    }
    
    //reset acceptances after equilibration
    for(auto &i : moves)
        i.reset_acceptance_counters();
    
    //estimators and vectors for data outputs
    Estimator estimator(local_comm);
    std::vector<std::vector<double> > energies;
    std::vector<std::vector<int> > windings;
    std::vector<std::vector<int> > permutations;
    std::vector<double> particles;
    energies.reserve(params.end);
    windings.reserve(params.end);
    permutations.reserve(params.end);
    particles.reserve(params.end);
    
    int counter = 0; //counts number of written data points
    bool path_dump = true; //boolean for whether to write out paths or not
    if(id == 0)
        std::cout << "Collecting Data ... " << std::endl;
    for(int step = 0; step < params.end; ++step){
        if(id == 0 && (step+1)%100 == 0)
            std::cout << "Simulation step " << step+1 << "/" <<params.end << ":\t" << params.real_particles << " particles" << std::endl;
        for(int j = 0; j < std::max(params.particles,params.init_particles); ++j){
            for(int n = 0; n < std::max(params.particles,params.init_particles); ++n){ //n^2 moves per step, once again
                int success = 0;
                do{ //one successful move per particle
                    int choice = 0;
                    if(params.worm_on) //choose appropriate move randomly
                        choice = move_choices[1][(rng.randint(move_choices[1].size()))];
                    else
                        choice = move_choices[0][(rng.randint(move_choices[0].size()))];
                    MPI_Bcast(&choice, 1, MPI_INT, 0, local_comm); //tell all processors about move
                    success = moves[choice].attempt(id, params, paths, rng, cos); //attempt move
                }while(success != 0);
            }
        }
        params.real_particles = params.particles - (ceil(double(params.worm_length)/params.total_slices) - double(params.worm_length)/params.total_slices);
        
        //if the configuration is diagonal, estimate thermodynamic quantities and write to file
        if(!params.worm_on && params.particles > 0){
            ++counter;
            if(counter%50 == 0 && path_dump) //path dump at certain intervals (can be changed)
                writer.path_dump(id, step, paths, params);
            std::vector<double> energy(6,0);
            std::vector<int> winding(params.dimensions,0);
            std::vector<int> permutation(params.particles,0);
            estimator.estimate(id, paths, params, energy, winding, permutation);
            energies.push_back(energy);
            windings.push_back(winding);
            permutations.push_back(permutation);
            if(id == 0)
                writer.write_step_state(step, params.particles, energy, permutation, winding);
        }
    }
    if(id == 0){
        std::cout << "Finalizing ... " << std::endl;
        writer.write_acceptances(counter, moves); //write final move acceptance ratios to output file
    }
}
