//
//  estimators.cpp
//  PIMC
//
//  Created by Adith Ramamurti on 10/25/16.
//  Copyright © 2016 Adith Ramamurti. All rights reserved.
//
#include "estimators.hpp"


/*---------------------------------------------------------------------------------------------------*

This file contains the implementation of the energy estimator for the PIMC simulation.

See estimators.hpp for class structure.

*---------------------------------------------------------------------------------------------------*/

//harmonic oscillator potential well (centered at zero)
inline double potential_value_harmonic(double dist2){
    double val = 0.5*dist2;
    return val;
}

//piecewise portion of Aziz potential
inline double aziz_pcws(double dist){
    if(dist >= 3.68335)
        return 1;
    else
        return exp(-pow((3.68335/dist-1),2));
}

//Aziz potential
inline double potential_value_aziz(double dist){
    double val = 10.8*(544850.4 * exp(-4.50018*dist)-(9424.94/pow(dist,10)+2556.63/pow(dist,8)+937.38/pow(dist,6))*aziz_pcws(dist));
    return val;
}

//Derivative of Aziz potential
inline double deriv_aziz_pcws(double dist){
    if(dist >= 3.69912)
        return 0;
    else
        return (7.3667*exp(-pow(3.68335/dist-1,2))*(3.68335/dist-1))/pow(dist,2);
}

//Gradient of Aziz potential
inline double grad_potential_value_aziz(double dist){
    double part1 = -2.45192E6*exp(-4.50018*dist);
    double part2 = (-94249.4/pow(dist,11)-20453./pow(dist,9)-5624.28/pow(dist,7))*aziz_pcws(dist);
    double part3 = (9424.94/pow(dist,10)+2556.63/pow(dist,8)+937.38/pow(dist,6))*deriv_aziz_pcws(dist);
    double val = 10.8*(part1-part2-part3);
    return val;
}

//Coulomb potential (computed using Ewald summation)
inline double potential_value_coulomb(std::vector<double>& dist, int chgi, int chgj, Parameters& params){
    if(params.coupling == 0) return 0;
    double val = 0;
    for(int nx = -params.nmax; nx <= params.nmax; nx++)
        for(int ny = -params.nmax; ny <= params.nmax; ny++)
            for(int nz = -params.nmax; nz <= params.nmax; nz++){
                double rSq = pow(dist[0] + nx*params.box_size, 2) + pow(dist[1] + ny*params.box_size, 2) + pow(dist[2] + nz*params.box_size, 2);
                double r = sqrt(rSq);
                if(rSq > 10E-10 && rSq < params.coulcut2){
                    val += chgi*chgj * erfc(params.alpha*r) / r;
                }
            }
    for(int kx = 0; kx <= params.kmax; kx++)
        for(int ky = 0; ky <= params.kmax; ky++)
            for(int kz = 0; kz <= params.kmax; kz++){
                std::vector<int> kVec = {kx, ky, kz};
                double k2 = 0;
                for(int dim = 0; dim < params.dimensions; ++dim)
                    k2 += pow(params.kfac * kVec[dim],2);
                if(k2 < params.kcut2 && (k2  != 0)){
                    int zeros = 0;
                    if(kx == 0)
                        ++zeros;
                    if(ky == 0)
                        ++zeros;
                    if(kz == 0)
                        ++zeros;
                    double sfac = 2.0;
                    if(zeros == 0)
                        sfac = 8.0;
                    else if(zeros == 1)
                        sfac = 4.0;
                    double efac = 1.0;
                    double expfac = 0;
                    for(int dim = 0; dim < params.dimensions; ++dim){
                        expfac = params.kfac * kVec[dim] * dist[dim];
                        efac = efac * cos(expfac);
                    }
                    double reci = sfac*(4*M_PI)/pow(params.box_size,3)*chgi*chgj*efac*exp(-k2*1/(4*params.alpha*params.alpha))/k2;
                    val += reci;
                }
            }
    if(std::abs(dist[0]) < 10E-10 && std::abs(dist[1]) < 10E-10 && std::abs(dist[2]) < 10E-10)
        val -= 2*params.alpha/sqrt(M_PI) * pow(chgi,2);
    val = params.coupling*val;
    return val;
}

//Gradient of Coulomb potential
inline std::vector<double> grad_potential_value_coulomb(std::vector<double>& dist, int chgi, int chgj, Parameters& params){
    if(params.coupling == 0) return std::vector<double>(params.dimensions, 0);
    std::vector<double> val(params.dimensions, 0);
    double gradV1 = 0;
    for(int nx = -params.nmax; nx <= params.nmax; nx++)
        for(int ny = -params.nmax; ny <= params.nmax; ny++)
            for(int nz = -params.nmax; nz <= params.nmax; nz++){
                double rSq = pow(dist[0] + nx*params.box_size, 2) + pow(dist[1] + ny*params.box_size, 2) + pow(dist[2] + nz*params.box_size, 2);
                double r = sqrt(rSq);
                if(rSq > 10E-10 && rSq < params.coulcut2){
                    gradV1 = chgi*chgj * (-2*exp(-rSq*pow(params.alpha,2))*params.alpha/(sqrt(M_PI)*r)- erfc(params.alpha*r) / rSq);
                    for(int dim = 0; dim < params.dimensions; ++dim){
                        val[dim] += gradV1 * dist[dim]/r;
                    }
                }
            }
    double gradV2 = 0;
    for(int kx = 0; kx <= params.kmax; kx++)
        for(int ky = 0; ky <= params.kmax; ky++)
            for(int kz = 0; kz <= params.kmax; kz++){
                std::vector<int> kVec = {kx, ky, kz};
                double k2 = 0;
                for(int dim = 0; dim < params.dimensions; ++dim)
                    k2 += pow(params.kfac * kVec[dim],2);
                if(k2 < params.kcut2 && (k2  != 0)){
                    int zeros = 0;
                    if(kx == 0)
                        ++zeros;
                    if(ky == 0)
                        ++zeros;
                    if(kz == 0)
                        ++zeros;
                    double sfac = 2.0;
                    if(zeros == 0)
                        sfac = 8.0;
                    else if(zeros == 1)
                        sfac = 4.0;
                    double efac = 1.0;
                    double expfac = 0;
                    for(int dim = 0; dim < params.dimensions; ++dim){
                        expfac = params.kfac * kVec[dim] * dist[dim];
                        efac = efac * sin(expfac);
                    }
                    gradV2 = -sfac*(4*M_PI)/pow(params.box_size,3)*chgi*chgj*efac*exp(-k2*1/(4*params.alpha*params.alpha))/k2;
                    for(int dim = 0; dim < params.dimensions; ++dim){
                        val[dim] += gradV2 * kVec[dim];
                    }
                }
            }
    for(int dim = 0; dim < params.dimensions; ++dim)
        val[dim] = params.coupling*val[dim];
    return val;
}


//Calculates energy, winding, and permutations for a configuration
void Estimator::estimate(int& id, Paths& paths, Parameters &params, std::vector<double>& energy, std::vector<int>& winding, std::vector<int>& permutations){
    std::vector<std::vector<std::vector<double> > > all_locations;
    
    //Gather bead locations from all processors (sent to master)
    if(id == 0){
        std::vector<std::vector<double> > slice;
        for(int i = 0; i < params.slices_per_process; ++i){
            for(int k = 0; k < params.particles; ++k){
                slice.push_back(paths.get_coordinate(i,k));
            }
            all_locations.push_back(slice);
            slice.clear();
        }
        std::vector<double> loc(params.dimensions);
        for(int j = 1; j < params.num_workers; ++j){
            for(int i = 0; i < params.slices_per_process; ++i){
                for(int k = 0; k < params.particles; ++k){
                    MPI_Recv(&loc[0], params.dimensions, MPI_DOUBLE, j, i*params.particles+k, local_comm, MPI_STATUS_IGNORE);
                    slice.push_back(loc);
                }
                all_locations.push_back(slice);
                slice.clear();
            }
        }
    }
    else{
        for(int i = 0; i < params.slices_per_process; ++i)
            for(int ptcl = 0; ptcl < params.particles; ++ptcl)
                MPI_Send(&paths.get_coordinate(i,ptcl)[0], params.dimensions, MPI_DOUBLE, 0, i*params.particles+ptcl, local_comm);
    }
    
    //Master calculates all values
    if(id == 0){
        double pe = 0; //potential energy
        double keth = 0; //kinetic (thermal) energy
        double kev = 0; // kinetic (virial) energy
        double kd = 0; //square kinetic distance
        
        //Set Coulomb potential Ewald summation parameters
        if(params.potential == 2){
            params.alpha = sqrt(M_PI)*pow(params.particles/params.volume2,1/6.);
            params.coulcut = sqrt(params.p)/params.alpha;
            params.coulcut2 = pow(params.coulcut,2);
            params.kcut = 2.*params.p/params.coulcut;
            params.kcut2 = pow(params.kcut,2);
            params.nmax = floor(params.coulcut/params.box_size);
            params.kmax = ceil(params.kcut/(2.*M_PI/params.box_size));
        }
        
        
        std::vector<double> dist(params.dimensions, 0);
        
        for(int s = 0; s < params.total_slices; ++s)
            for(int p1 = 0; p1 < params.particles; ++p1){
                //For each slice, calculate potential
                switch(params.potential){
                    case 0:
                        pe += potential_value_harmonic(inner_product(all_locations[s][p1].begin(), all_locations[s][p1].end(), all_locations[s][p1].begin(), 0.0));
                        break;
                    case 1:
                        for (int p2 = p1+1; p2 < params.particles; ++p2){
                            distance(all_locations[s][p2] , all_locations[s][p1],dist, params.box_size);
                            double pd = sqrt(inner_product(dist.begin(),dist.end(),dist.begin(),0.0));
                            pe += potential_value_aziz(pd);
                        }
                        break;
                    case 2:
                        for (int p2 = 0; p2 < params.particles; ++p2){
                            distance(all_locations[s][p2] , all_locations[s][p1],dist, params.box_size);
                            pe += potential_value_coulomb(dist, paths.charge[p2], paths.charge[p1], params);
                        }
                        break;
                }
                
                //and the square distance between each bead
                if(s+1 != params.total_slices){
                    distance(all_locations[s+1][p1] , all_locations[s][p1],dist, params.box_size);
                    kd += inner_product(dist.begin(),dist.end(),dist.begin(),0.0);
                }
                else{
                    distance(all_locations[0][paths.forward_connects[p1]] , all_locations[s][p1],dist, params.box_size);
                    kd += inner_product(dist.begin(),dist.end(),dist.begin(),0.0);
                }
            }
        if(params.potential == 2)
            pe = (pe/2);
        pe = pe/params.total_slices;
        
        //Kinetic (thermal) energy
        keth = 0.5*params.dimensions*params.particles/params.tau  - 1/(4.*params.lambda*pow(params.tau,2))*kd/params.total_slices;
        
        std::vector<std::vector<std::vector<double> > > trajectories;
        std::vector<double> total_distance(params.dimensions, 0.0);
        trajectories.reserve(all_locations.size());
        trajectories.push_back(all_locations[0]);
        std::vector<std::vector<double> > slice;
        
        //calculate total distances and positions for each bead on each path (accounting for periodic conditions)
        for(int j = 1; j < params.total_slices; ++j){
            for(int i = 0; i < params.particles; ++i){
                distance(all_locations[j-1][i], all_locations[j][i], dist, params.box_size);
                std::vector<double> new_loc(params.dimensions,0);
                for(int d = 0; d < params.dimensions; ++d){
                    new_loc[d] = trajectories.back()[i][d]+dist[d];
                    total_distance[d] += dist[d];
                }
                slice.push_back(new_loc);
            }
            trajectories.push_back(slice);
            slice.clear();
        }
        for(int i = 0; i < params.particles; ++i){
            distance(all_locations[params.total_slices-1][i], all_locations[0][paths.forward_connects[i]], dist, params.box_size);
            std::vector<double> new_loc(params.dimensions,0);
            for(int d = 0; d < params.dimensions; ++d){
                new_loc[d] = trajectories.back()[i][d]+dist[d];
                total_distance[d] += dist[d];
            }
            slice.push_back(new_loc);
        }
        trajectories.push_back(slice);
        slice.clear();
        for(int j = 1; j < params.total_slices; ++j){
            for(int i = 0; i < params.particles; ++i){
                distance(all_locations[j-1][paths.forward_connects[i]], all_locations[j][paths.forward_connects[i]], dist, params.box_size);
                std::vector<double> new_loc(params.dimensions,0);
                for(int d = 0; d < params.dimensions; ++d)
                    new_loc[d] = trajectories.back()[i][d]+dist[d];
                slice.push_back(new_loc);
            }
            trajectories.push_back(slice);
            slice.clear();
        }
        
        //Calculate energy using virial method
        
        //set constants
        int virial_window = params.total_slices;
        double E1 = 0.5*params.dimensions*params.particles/(params.tau*virial_window);
        double E2fac = 1/(4*params.lambda*params.total_slices*virial_window*params.tau*params.tau);
        double E2 = 0;
        double E3fac = 1./(2*params.total_slices);
        double E3 = 0;
        
        //virial window distances
        for(int ptcl = 0; ptcl < params.particles; ++ptcl){
            for(int slice = 0; slice < params.total_slices; ++slice){
                std::vector<double> v1;
                std::vector<double> v2;
                v1.reserve(params.dimensions);
                std::transform(trajectories[slice+virial_window][ptcl].begin(), trajectories[slice+virial_window][ptcl].end(), trajectories[slice][ptcl].begin(),
                               std::back_inserter(v1), std::minus<double>());
                v2.reserve(params.dimensions);
                std::transform(trajectories[slice+virial_window-1][ptcl].begin(), trajectories[slice+virial_window-1][ptcl].end(), trajectories[slice+virial_window][ptcl].begin(), std::back_inserter(v2), std::minus<double>());
                E2 += inner_product(v1.begin(), v1.end(), v2.begin(), 0.0);
            }
        }
        
        //center of mass offsets
        std::vector<std::vector<double> > com(params.particles);
        for(int ptcl = 0; ptcl < params.particles; ++ptcl){
            com[ptcl].resize(params.dimensions,0);
            for(int slice = 0; slice < params.total_slices; ++slice){
                for(int dim = 0; dim < params.dimensions; ++dim){
                    com[ptcl][dim] += trajectories[slice][ptcl][dim];
                }
            }
            for(int dim = 0; dim < params.dimensions; ++dim)
                com[ptcl][dim] = com[ptcl][dim]/params.total_slices;
        }
        for(int p1 = 0; p1 < params.particles; ++p1)
            for(int s = 0; s < params.total_slices; s++)
                switch(params.potential){
                    case 0:
                        for(int dim = 0; dim < params.dimensions; ++dim){
                            E3 += (all_locations[p1][s][dim] - com[p1][dim])*all_locations[p1][s][dim];
                        }
                        break;
                    case 1:
                        for(int p2 = 0; p2 < params.particles; ++p2)
                            if(p1 != p2){
                                distance(all_locations[s][p2],all_locations[s][p1],dist, params.box_size);
                                double npd = sqrt(inner_product(dist.begin(),dist.end(),dist.begin(),0.0));
                                std::vector<double> pd = dist;
                                std::vector<double> gradVvec(params.dimensions, 0);
                                double gradV = grad_potential_value_aziz(npd);
                                for(int d = 0; d < params.dimensions; ++d){
                                    gradVvec[d] = gradV*pd[d]/npd;
                                }
                                distance(com[p1],all_locations[s][p1],dist, params.box_size);
                                E3 += inner_product(dist.begin(), dist.end(), gradVvec.begin(), 0.0);
                            }
                        break;
                    case 2:
                        for(int p2 = 0; p2 < params.particles; ++p2){
                            distance(all_locations[s][p2],all_locations[s][p1],dist, params.box_size);
                            std::vector<double> gradVvec = grad_potential_value_coulomb(dist, paths.charge[p2], paths.charge[p1], params);
                            distance(com[p1],all_locations[s][p1],dist, params.box_size);
                            E3 += inner_product(dist.begin(), dist.end(), gradVvec.begin(), 0.0);
                        }
                        break;
                }
        //total virial energy
        kev = E1 + E2fac * E2 + E3fac * E3;
        
        //calculate permutations based on forward connects values
        for(int part = 0; part < params.particles; ++part){
            int counter = 0;
            int particle = paths.forward_connects[part];
            while(particle != part){
                ++counter;
                particle = paths.forward_connects[particle];
            }
            ++permutations[counter];
        }
        
        //calculate winding
        for(int dim = 0; dim < params.dimensions; ++dim)
            winding[dim] = round(total_distance[dim]/params.box_size);
        
        //combine energy calculations into a vector
        energy = {pe+keth,keth,pe, kev+pe, kev, pe};        
    }
}
