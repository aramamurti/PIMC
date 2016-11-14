//
//  moves.cpp
//  PIMC
//
//  Created by Adith Ramamurti on 10/17/16.
//  Copyright © 2016 Adith Ramamurti. All rights reserved.
//

#include <stdio.h>
#include "moves.hpp"

inline double potential_value_harmonic(double dist2){
    double val = 0.5*dist2;
    return val;
}

inline double aziz_pcws(double dist){
    if(dist >= 3.68335)
        return 1;
    else
        return exp(-pow((3.68335/dist-1),2));
}

inline double potential_value_aziz(double dist){
    double val = 10.8*(544850.4 * exp(-4.50018*dist)-(9424.94/pow(dist,10)+2556.63/pow(dist,8)+937.38/pow(dist,6))*aziz_pcws(dist));
    return val;
}

inline double potential_value_coulomb(std::vector<double>& dist, std::vector<int>& kVec, int chgi, int chgj, Parameters& params){
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
                kVec = {kx, ky, kz};
                double k2 = 0;
                for(int dim = 0; dim < params.dimensions; ++dim)
                    k2 += pow(params.kfac * kVec[dim],2);
                if(k2 < params.kcut2 && (k2  != 0)){
                    int zeros = 0;
                    if(kx == 0)
                        zeros++;
                    if(ky == 0)
                        zeros++;
                    if(kz == 0)
                        zeros++;
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
    return params.coupling*val;
}


Moves::Moves(MPI_Comm &local){
    num_attempts = 0;
    num_accepts = 0;
    local_comm = local;
}

int Moves::attempt(int &id, Parameters &params, Paths &paths, RNG &rng){
    return 0;
}

void Moves::reset_acceptance_counters(){
    num_attempts = 0;
    num_accepts = 0;
}

Center_of_Mass::Center_of_Mass(int &id, Parameters &params, MPI_Comm &local) : Moves(local){
    delta = sqrt(params.lambda*params.tau);
    move_name = "Center of Mass";
    worm_off = true;
    worm_on = true;
}

int Center_of_Mass::attempt(int &id,Parameters &params, Paths &paths, RNG &rng){
    Moves::attempt(id, params, paths, rng);
    if(params.particles == 0) return 0;
    std::vector<double> shift(params.dimensions);
    for(auto &dim : shift)
        dim = rng.randgaussian(delta);
    int ptcl = rng.randint(params.particles);
    MPI_Bcast(&shift[0], params.dimensions, MPI_DOUBLE, 0, local_comm);
    MPI_Bcast(&ptcl, 1, MPI_INT, 0, local_comm);
    int cur_part = ptcl;
    int end_part = ptcl;
    if(paths.broken[ptcl]){
        cur_part = params.worm_head.second;
        end_part = -1;
    }
    do{
        ptcls.push_back(cur_part);
        cur_part = paths.forward_connects[cur_part];
    }while(cur_part != end_part);
    new_coordinates.resize(params.slices_per_process*ptcls.size());
    new_coordinates_ahead.resize(new_coordinates.size());
    keys_ahead.resize(new_coordinates.size());
    new_distances.reserve(ptcls.size()*params.particles*params.slices_per_process);
    for(int i = 0; i < ptcls.size(); ++i){
        for(int slice = 0; slice < params.slices_per_process; ++slice){
            if(!(ptcls[i] == params.worm_head.second && slice+params.my_start < params.worm_head.first) && !(ptcls[i] == params.worm_tail.second && slice+params.my_start > params.worm_tail.first))
                std::transform(paths.get_coordinate(slice,ptcls[i]).begin(),paths.get_coordinate(slice,ptcls[i]).end(),shift.begin(),std::back_inserter(new_coordinates[i*params.slices_per_process+slice]), std::plus<double>());
        }
    }
    check(id, params, paths, rng);
    MPI_Bcast(&ac_re, 1, MPI_INT, 0, local_comm);
    num_attempts++;
    if(ac_re){
        num_accepts++;
        for(int i = 0; i < ptcls.size(); ++i)
            for(int slice = 0; slice < params.slices_per_process; ++slice)
                if(!(ptcls[i] == params.worm_head.second && slice+params.my_start < params.worm_head.first) && !(ptcls[i] == params.worm_tail.second && slice+params.my_start > params.worm_tail.first)){
                    put_in_box(new_coordinates[i*params.slices_per_process+slice], params.box_size);
                    paths.set_coordinate(slice, ptcls[i],new_coordinates[i*params.slices_per_process+slice], false);
                }
        paths.update_separations(new_distances);
        if(!params.gce){
            for(int i = 0; i < ptcls.size(); ++i)
                for(int slice = 0; slice < params.slices_per_process; ++slice)
                    if(!(ptcls[i] == params.worm_head.second && slice+params.my_start < params.worm_head.first) && !(ptcls[i] == params.worm_tail.second && slice+params.my_start > params.worm_tail.first))
                        paths.calculate_kinetic_separations_start(slice, paths.get_coordinate_key(slice, ptcls[i]));
            for(int i = 0; i < ptcls.size(); ++i)
                for(int tslice = 0; tslice < params.total_slices; ++tslice){
                    int slice_back = positive_modulo(tslice-params.multistep_dist, params.total_slices);
                    if(tslice >= params.my_start && tslice <= params.my_end){
                        int slice = tslice%params.slices_per_process;
                        int send_rank = slice_back/params.slices_per_process;
                        int send_tag = i*params.total_slices+slice_back;
                        new_coordinates[i*params.slices_per_process+slice].resize(params.dimensions);
                        int key = paths.get_coordinate_key(slice, ptcls[i]);
                        if(send_rank != id){
                            MPI_Send(&new_coordinates[i*params.slices_per_process+slice][0], params.dimensions, MPI_DOUBLE, send_rank, send_tag, local_comm);
                            MPI_Send(&key, 1, MPI_INT, send_rank, params.total_slices*ptcls.size()+send_tag, local_comm);
                        }
                        else{
                            keys_ahead[i*params.slices_per_process+slice_back%params.slices_per_process] = key;
                            new_coordinates_ahead[i*params.slices_per_process+slice_back%params.slices_per_process] = new_coordinates[i*params.slices_per_process+slice];
                            paths.set_kinetic_end(slice_back%params.slices_per_process, keys_ahead[i*params.slices_per_process+slice_back%params.slices_per_process], new_coordinates_ahead[i*params.slices_per_process+slice_back%params.slices_per_process], true);
                        }
                    }
                    else if(slice_back >= params.my_start && slice_back <= params.my_end){
                        int slice = slice_back%params.slices_per_process;
                        int recv_rank = tslice/params.slices_per_process;
                        int recv_tag = i*params.total_slices+slice_back;
                        new_coordinates_ahead[i*params.slices_per_process+slice].resize(params.dimensions);
                        MPI_Recv(&new_coordinates_ahead[i*params.slices_per_process+slice][0], params.dimensions, MPI_DOUBLE, recv_rank, recv_tag, local_comm, MPI_STATUS_IGNORE);
                        MPI_Recv(&keys_ahead[i*params.slices_per_process+slice], 1, MPI_INT, recv_rank, params.total_slices*ptcls.size()+recv_tag, local_comm, MPI_STATUS_IGNORE);
                        if(keys_ahead[i*params.slices_per_process+slice] != 0){
                            paths.set_kinetic_end(slice, keys_ahead[i*params.slices_per_process+slice], new_coordinates_ahead[i*params.slices_per_process+slice], true);
                        }
                    }
                }
        }
    }
    ptcls.clear();
    new_coordinates.clear();
    new_coordinates_ahead.clear();
    keys_ahead.clear();
    new_distances.clear();
    return 0;
}

void Center_of_Mass::check(int &id, Parameters &params, Paths &paths, RNG &rng){
    std::vector<double> old_action_p(params.slices_per_process,0);
    std::vector<double> new_action_p(params.slices_per_process,0);
    std::vector<double> dist(params.dimensions);
    std::vector<int> kVec(params.dimensions);
    if(params.potential == 2){
        params.alpha = sqrt(M_PI)*pow(params.particles/params.volume2,1/6.);
        params.coulcut = sqrt(params.p)/params.alpha;
        params.coulcut2 = pow(params.coulcut,2);
        params.kcut = 2.*params.p/params.coulcut;
        params.kcut2 = pow(params.kcut,2);
        params.nmax = floor(params.coulcut/params.box_size);
        params.kmax = ceil(params.kcut/(2.*M_PI/params.box_size));
    }
    for(int j = 0; j < ptcls.size(); ++j){
        for(int i = 0; i < params.slices_per_process; ++i){
            if(!(ptcls[j] == params.worm_head.second && i+params.my_start < params.worm_head.first) && !(ptcls[j] == params.worm_tail.second && i+params.my_start > params.worm_tail.first)){
                double old_action_local = 0;
                double new_action_local = 0;
                switch(params.potential){
                    case 0:
                        new_action_p[i]+= potential_value_harmonic(inner_product(new_coordinates[j*params.slices_per_process+i].begin(),new_coordinates[j*params.slices_per_process+i].end(),new_coordinates[j*params.slices_per_process+i].begin(),0.0));
                        old_action_p[i] = potential_value_harmonic(inner_product(paths.get_coordinate(i,ptcls[j]).begin(),paths.get_coordinate(i,ptcls[j]).end(),paths.get_coordinate(i,ptcls[j]).begin(),0.0));
                        break;
                    case 1:
#pragma omp parallel for firstprivate(dist) reduction (+:new_action_local, old_action_local)
                        for(int ptcl = 0; ptcl < params.particles; ++ptcl){
                            if(!(ptcl == params.worm_head.second && i+params.my_start < params.worm_head.first) && !(ptcl == params.worm_tail.second && i+params.my_start > params.worm_tail.first)){
                                if(std::find(ptcls.begin(), ptcls.end(), ptcl) == ptcls.end()){
                                    old_action_local += potential_value_aziz(paths.get_separation(i, ptcl, ptcls[j]));
                                    distance(paths.get_coordinate(i,ptcl), new_coordinates[j*params.slices_per_process+i], dist, params.box_size);
                                    double r = sqrt(std::inner_product(dist.begin(), dist.end(), dist.begin(), 0.0));
                                    new_action_local += potential_value_aziz(r);
#pragma omp critical
                                    new_distances.push_back(std::tuple<std::pair<int,int>,std::vector<double>,double>(std::pair<int,int>(paths.get_coordinate_key(i, ptcl), paths.get_coordinate_key(i, ptcls[j])), dist, r));
                                }
                            }
                        }
                        old_action_p[i] += old_action_local;
                        new_action_p[i] += new_action_local;
                        break;
                    case 2:
#pragma omp parallel for firstprivate(dist) reduction (+:new_action_local, old_action_local)
                        for(int ptcl = 0; ptcl < params.particles; ++ptcl){
                            if(!(ptcl == params.worm_head.second && i+params.my_start < params.worm_head.first) && !(ptcl == params.worm_tail.second && i+params.my_start > params.worm_tail.first)){
                                if(std::find(ptcls.begin(), ptcls.end(), ptcl) == ptcls.end()){
                                    old_action_local += potential_value_coulomb(paths.get_separation_vector(i, ptcl, ptcls[j]), kVec, paths.charge[ptcl],paths.charge[ptcls[j]],  params);
                                    double r = 0;
                                    distance(paths.get_coordinate(i,ptcl), new_coordinates[j*params.slices_per_process+i], dist, params.box_size);
                                    r = sqrt(std::inner_product(dist.begin(), dist.end(), dist.begin(), 0.0));
                                    new_action_local += potential_value_coulomb(dist, kVec,  paths.charge[ptcl],paths.charge[ptcls[j]], params);
#pragma omp critical
                                    new_distances.push_back(std::tuple<std::pair<int,int>,std::vector<double>,double>(std::pair<int,int>(paths.get_coordinate_key(i, ptcl), paths.get_coordinate_key(i, ptcls[j])), dist, r));
                                }
                            }
                        }
                        old_action_p[i] += old_action_local;
                        new_action_p[i] += new_action_local;
                        break;
                }
            }
        }
    }
    if(id != 0){
        MPI_Send(&old_action_p[0], params.slices_per_process, MPI_DOUBLE, 0, 0, local_comm);
        MPI_Send(&new_action_p[0], params.slices_per_process, MPI_DOUBLE, 0, 1, local_comm);
    }
    else{
        ac_re = false;
        std::vector<double> old_action;
        std::vector<double> new_action;
        old_action.reserve(params.total_slices);
        new_action.reserve(params.total_slices);
        old_action.insert(old_action.end(),old_action_p.begin(),old_action_p.end());
        new_action.insert(new_action.end(),new_action_p.begin(),new_action_p.end());
        for(int j = 1; j < params.num_workers; ++j){
            MPI_Recv(&old_action_p[0], params.slices_per_process, MPI_DOUBLE, j, 0, local_comm, MPI_STATUS_IGNORE);
            MPI_Recv(&new_action_p[0], params.slices_per_process, MPI_DOUBLE, j, 1, local_comm, MPI_STATUS_IGNORE);
            old_action.insert(old_action.end(),old_action_p.begin(),old_action_p.end());
            new_action.insert(new_action.end(),new_action_p.begin(),new_action_p.end());
        }
        double oa = 0;
        double na = 0;
        for(int j = 0; j < params.total_slices; ++j){
            oa += params.tau*(old_action[j]);
            na += params.tau*(new_action[j]);
        }
        if(rng.randnormed(1) < exp(-(na - oa)))
            ac_re = true;
    }
}

Bisection::Bisection(int& id, Parameters &params, MPI_Comm &local) : Moves(local){
    multistep_dist = params.multistep_dist;
    minp.resize(params.num_workers);
    for(int i = 0; i < params.total_slices; ++i){
        multisteps.push_back(std::vector<int>(1,i));
        minp[i/params.slices_per_process].push_back(i);
        for(int j = 1; j <= multistep_dist; ++j){
            multisteps[i].push_back((i+j)%params.total_slices);
            minp[((i+j)%params.total_slices)/params.slices_per_process].push_back(i);
        }
    }
    for(auto &i : minp){
        sort(i.begin(), i.end());
        auto last = unique(i.begin(), i.end());
        i.erase(last, i.end());
    }
    move_name = "Bisection";
    worm_off = true;
    worm_on = true;
}

int Bisection::attempt(int &id,Parameters &params, Paths &paths, RNG &rng){
    Moves::attempt(id, params, paths, rng);
    if(params.particles == 0) return 0;
    int set = rng.randint(multisteps.size());
    bisection_info.clear();
    bisection_info.push_back(multisteps[set].front());
    bisection_info.push_back(multisteps[set].back());
    if(multisteps[set].back() < multisteps[set].front())
        bisection_info.back() += params.total_slices;
    bisection_info.push_back(rng.randint(params.particles));
    MPI_Bcast(&bisection_info[0], 3, MPI_INT, 0,local_comm);
    if(paths.broken[bisection_info[2]]){
        int final_ptcl = bisection_info[2];
        if(bisection_info[1] > params.total_slices)
            final_ptcl = paths.forward_connects[final_ptcl];
        if((bisection_info[2] == params.worm_head.second && bisection_info[0] < params.worm_head.first) || (bisection_info[2] == params.worm_tail.second && bisection_info[1] > params.worm_tail.first)|| (final_ptcl == params.worm_tail.second && bisection_info[0]+multistep_dist > params.worm_tail.first)){
            return 1;
        }
    }
    new_coordinates.resize(params.slices_per_process);
    new_coordinates_ahead.resize(params.slices_per_process);
    keys_ahead.resize(params.slices_per_process);
    new_distances.reserve(params.slices_per_process*params.particles);
    int level = int(round(log2(bisection_info[1]-bisection_info[0])));
    int counter = 0;
    ptcl_slice.resize(params.slices_per_process,bisection_info[2]);
    if(bisection_info[1]>= params.total_slices)
        for(int i = bisection_info[0]; i <=bisection_info[1]; ++i)
            if(i%params.total_slices >= params.my_start && i%params.total_slices <= params.my_end && i >= params.total_slices)
                ptcl_slice[i%params.slices_per_process] = paths.forward_connects[bisection_info[2]];
    while(level){
        int num_comps = pow(2,counter);
        counter++;
        std::vector<int> slices_to_compute;
        for(int i = 0; i < num_comps; ++i)
            slices_to_compute.push_back(int((1/pow(2.,counter) + i/pow(2.,counter-1))*(bisection_info[1]-bisection_info[0])+bisection_info[0])%params.total_slices);
        for(auto &slice : slices_to_compute){
            int first = int(slice - pow(2,level-1)+params.total_slices)%params.total_slices;
            int second = int(slice + pow(2, level-1))%params.total_slices;
            if(slice >= params.my_start && slice <= params.my_end){
                std::vector<double> start(params.dimensions);
                std::vector<double> end(params.dimensions);
                if(first < params.my_start || first > params.my_end)
                    MPI_Recv(&start[0], params.dimensions,MPI_DOUBLE,first/params.slices_per_process,2*counter-1,local_comm,MPI_STATUS_IGNORE);
                else{
                    if(counter == 1){
                        start = paths.get_coordinate(first%params.slices_per_process, ptcl_slice[first%params.slices_per_process]);
                        new_coordinates[first%params.slices_per_process] = start;
                    }
                    else
                        start = new_coordinates[first%params.slices_per_process];
                }
                if(second < params.my_start || second > params.my_end)
                    MPI_Recv(&end[0], params.dimensions,MPI_DOUBLE,second/params.slices_per_process,2*counter,local_comm,MPI_STATUS_IGNORE);
                else{
                    if(counter == 1){
                        end = paths.get_coordinate(second%params.slices_per_process, ptcl_slice[second%params.slices_per_process]);
                        new_coordinates[second%params.slices_per_process] = end;
                    }
                    else
                        end = new_coordinates[second%params.slices_per_process];
                }
                average_loc(start, end, new_coordinates[slice%params.slices_per_process], params.box_size);
                double tau = params.tau*pow(2, level-1);
                double width = sqrt(params.lambda*tau);
                for(auto &dim : new_coordinates[slice%params.slices_per_process]){
                    dim += rng.randgaussian(width);
                }
            }
            else{
                if(first >= params.my_start && first <= params.my_end){
                    if(counter == 1)
                        new_coordinates[first%params.slices_per_process] = paths.get_coordinate(first%params.slices_per_process, ptcl_slice[first%params.slices_per_process]);
                    MPI_Send(&new_coordinates[first%params.slices_per_process][0],params.dimensions,MPI_DOUBLE,slice/params.slices_per_process,2*counter-1,local_comm);
                }
                if(second >= params.my_start && second <= params.my_end){
                    if(counter == 1)
                        new_coordinates[second%params.slices_per_process] = paths.get_coordinate(second%params.slices_per_process, ptcl_slice[second%params.slices_per_process]);
                    MPI_Send(&new_coordinates[second%params.slices_per_process][0],params.dimensions,MPI_DOUBLE,slice/params.slices_per_process,2*counter,local_comm);
                }
            }
        }
        --level;
    }
    check(id, params, paths, rng);
    MPI_Bcast(&ac_re, 1, MPI_INT, 0, local_comm);
    ++num_attempts;
    if(ac_re){
        ++num_accepts;
        for(int i = bisection_info[0]; i <= bisection_info[1]; ++i)
            if(i%params.total_slices >= params.my_start && i%params.total_slices <= params.my_end){
                int slice = (i%params.total_slices)%params.slices_per_process;
                put_in_box(new_coordinates[slice], params.box_size);
                paths.set_coordinate(slice, ptcl_slice[slice],new_coordinates[slice], false);
            }
        paths.update_separations(new_distances);
        if(!params.gce){
            for(int i = bisection_info[0]; i <= bisection_info[1]; ++i)
                if(i%params.total_slices >= params.my_start && i%params.total_slices <= params.my_end){
                    int slice = (i%params.total_slices)%params.slices_per_process;
                    paths.calculate_kinetic_separations_start(slice, paths.get_coordinate_key(slice, ptcl_slice[slice]));
                }
            for(int i = bisection_info[0]; i <= bisection_info[1]; ++i){
                int slice_back = positive_modulo(i-params.multistep_dist, params.total_slices);
                int tslice = i%params.total_slices;
                if(tslice >= params.my_start && tslice <= params.my_end){
                    int slice = tslice%params.slices_per_process;
                    int send_rank = slice_back/params.slices_per_process;
                    int send_tag = slice_back;
                    new_coordinates[slice].resize(params.dimensions);
                    int key = paths.get_coordinate_key(slice, ptcl_slice[slice]);
                    if(slice_back < params.my_start || slice_back > params.my_end){
                        MPI_Send(&new_coordinates[slice][0], params.dimensions, MPI_DOUBLE, send_rank, send_tag, local_comm);
                        MPI_Send(&key, 1, MPI_INT, send_rank, send_tag+params.slices_per_process, local_comm);
                    }
                    else{
                        keys_ahead[slice_back%params.slices_per_process] = key;
                        new_coordinates_ahead[slice_back%params.slices_per_process] = new_coordinates[slice];
                        paths.set_kinetic_end(slice_back%params.slices_per_process, keys_ahead[slice_back%params.slices_per_process], new_coordinates_ahead[slice_back%params.slices_per_process]);
                    }
                }
                else if(slice_back >= params.my_start && slice_back <= params.my_end){
                    int slice = slice_back%params.slices_per_process;
                    int recv_rank = tslice/params.slices_per_process;
                    int recv_tag = slice_back;
                    new_coordinates_ahead[slice].resize(params.dimensions);
                    MPI_Recv(&new_coordinates_ahead[slice][0], params.dimensions, MPI_DOUBLE, recv_rank, recv_tag, local_comm, MPI_STATUS_IGNORE);
                    MPI_Recv(&keys_ahead[slice], 1, MPI_INT, recv_rank, recv_tag+params.slices_per_process, local_comm, MPI_STATUS_IGNORE);
                    paths.set_kinetic_end(slice, keys_ahead[slice], new_coordinates_ahead[slice]);
                }
            }
        }
    }
    bisection_info.clear();
    ptcl_slice.clear();
    new_coordinates.clear();
    new_coordinates_ahead.clear();
    keys_ahead.clear();
    new_distances.clear();
    return 0;
}

void Bisection::check(int &id, Parameters &params, Paths &paths, RNG &rng){
    std::vector<double> old_action_p(params.slices_per_process,0);
    std::vector<double> new_action_p(params.slices_per_process,0);
    std::vector<double> dist(params.dimensions);
    std::vector<int> kVec(params.dimensions);
    if(params.potential == 2){
        params.alpha = sqrt(M_PI)*pow(params.particles/params.volume2,1/6.);
        params.coulcut = sqrt(params.p)/params.alpha;
        params.coulcut2 = pow(params.coulcut,2);
        params.kcut = 2.*params.p/params.coulcut;
        params.kcut2 = pow(params.kcut,2);
        params.nmax = floor(params.coulcut/params.box_size);
        params.kmax = ceil(params.kcut/(2.*M_PI/params.box_size));
    }
    for(int i = bisection_info[0]; i <= bisection_info[1]; ++i){
        int slice =i%params.total_slices;
        if(slice >= params.my_start && slice <= params.my_end){
            double old_action_local = 0;
            double new_action_local = 0;
            switch(params.potential){
                case 0:
                    old_action_p[slice%params.slices_per_process] += potential_value_harmonic(inner_product(paths.get_coordinate(slice%params.slices_per_process,ptcl_slice[slice%params.slices_per_process]).begin(),paths.get_coordinate(slice%params.slices_per_process,ptcl_slice[slice%params.slices_per_process]).end(),paths.get_coordinate(slice%params.slices_per_process,ptcl_slice[slice%params.slices_per_process]).begin(),0.0));
                    new_action_p[slice%params.slices_per_process] += potential_value_harmonic(inner_product(new_coordinates[slice%params.slices_per_process].begin(),new_coordinates[slice%params.slices_per_process].end(),new_coordinates[slice%params.slices_per_process].begin(),0.0));
                    break;
                case 1:
#pragma omp parallel for firstprivate(dist) reduction (+:new_action_local, old_action_local)
                    for(int ptcl = 0; ptcl < params.particles; ++ptcl){
                        if(!(ptcl == params.worm_head.second && slice < params.worm_head.first) && !(ptcl == params.worm_tail.second && slice > params.worm_tail.first)){
                            if(ptcl != ptcl_slice[slice%params.slices_per_process]){
                                old_action_local += potential_value_aziz(paths.get_separation(slice%params.slices_per_process, ptcl, ptcl_slice[slice%params.slices_per_process]));
                                distance(paths.get_coordinate(slice%params.slices_per_process,ptcl), new_coordinates[slice%params.slices_per_process], dist, params.box_size);
                                double r = sqrt(std::inner_product(dist.begin(), dist.end(), dist.begin(), 0.0));
                                new_action_local += potential_value_aziz(r);
#pragma omp critical
                                new_distances.push_back(std::tuple<std::pair<int,int>,std::vector<double>,double>(std::pair<int,int>(paths.get_coordinate_key(slice%params.slices_per_process, ptcl),paths.get_coordinate_key(slice%params.slices_per_process, ptcl_slice[slice%params.slices_per_process])), dist, r));
                            }
                        }
                    }
                    old_action_p[slice%params.slices_per_process] += old_action_local;
                    new_action_p[slice%params.slices_per_process] += new_action_local;
                    break;
                case 2:
#pragma omp parallel for firstprivate(dist) reduction (+:new_action_local, old_action_local)
                    for(int ptcl = 0; ptcl < params.particles; ++ptcl){
                        if(!(ptcl == params.worm_head.second && slice < params.worm_head.first) && !(ptcl == params.worm_tail.second && slice > params.worm_tail.first)){
                            if(ptcl != ptcl_slice[slice%params.slices_per_process]){
                                old_action_local += potential_value_coulomb(paths.get_separation_vector(slice%params.slices_per_process, ptcl, ptcl_slice[slice%params.slices_per_process]), kVec,  paths.charge[ptcl],paths.charge[ptcl_slice[slice%params.slices_per_process]], params);
                                distance(paths.get_coordinate(slice%params.slices_per_process,ptcl), new_coordinates[slice%params.slices_per_process], dist, params.box_size);
                                double r = sqrt(std::inner_product(dist.begin(), dist.end(), dist.begin(), 0.0));
                                new_action_local += potential_value_coulomb(dist, kVec, paths.charge[ptcl], paths.charge[ptcl_slice[slice%params.slices_per_process]],  params);
#pragma omp critical
                            new_distances.push_back(std::tuple<std::pair<int,int>,std::vector<double>,double>(std::pair<int,int>(paths.get_coordinate_key(slice%params.slices_per_process, ptcl),paths.get_coordinate_key(slice%params.slices_per_process, ptcl_slice[slice%params.slices_per_process])), dist, r));
                            }
                        }
                    }
                    old_action_p[slice%params.slices_per_process] = old_action_local;
                    new_action_p[slice%params.slices_per_process] = new_action_local;
                    break;
            }
        }
    }
    if(id != 0){
        MPI_Send(&old_action_p[0], params.slices_per_process, MPI_DOUBLE, 0, 0, local_comm);
        MPI_Send(&new_action_p[0], params.slices_per_process, MPI_DOUBLE, 0, 1, local_comm);
    }
    else{
        ac_re = false;
        std::vector<double> old_action;
        std::vector<double> new_action;
        old_action.reserve(params.total_slices);
        new_action.reserve(params.total_slices);
        old_action.insert(old_action.end(),old_action_p.begin(),old_action_p.end());
        new_action.insert(new_action.end(),new_action_p.begin(),new_action_p.end());
        for(int i = 1; i < params.num_workers; ++i){
            MPI_Recv(&old_action_p[0], params.slices_per_process, MPI_DOUBLE, i, 0, local_comm, MPI_STATUS_IGNORE);
            MPI_Recv(&new_action_p[0], params.slices_per_process, MPI_DOUBLE, i, 1, local_comm, MPI_STATUS_IGNORE);
            old_action.insert(old_action.end(),old_action_p.begin(),old_action_p.end());
            new_action.insert(new_action.end(),new_action_p.begin(),new_action_p.end());
        }
        double oa = 0;
        double na = 0;
        for(int j = bisection_info[0]; j < bisection_info[1]; ++j){
            int slice = j%params.total_slices;
            int slicep1 = (j+1)%params.total_slices;
            oa += params.tau/2.*(old_action[slice] + old_action[slicep1]);
            na += params.tau/2.*(new_action[slice] + new_action[slicep1]);
        }
        if(rng.randnormed(1) < exp(-(na - oa)))
            ac_re = true;
    }
}

Permutation_Bisection::Permutation_Bisection(int& id, Parameters &params, MPI_Comm &local) : Moves(local){
    multistep_dist = params.multistep_dist;
    minp.resize(params.num_workers);
    for(int i = 0; i < params.total_slices; ++i){
        multisteps.push_back(std::vector<int>(1,i));
        minp[i/params.slices_per_process].push_back(i);
        for(int j = 1; j <= multistep_dist; ++j){
            multisteps[i].push_back((i+j)%params.total_slices);
            minp[((i+j)%params.total_slices)/params.slices_per_process].push_back(i);
        }
    }
    for(auto &i : minp){
        sort(i.begin(), i.end());
        auto last = unique(i.begin(), i.end());
        i.erase(last, i.end());
    }
    move_name = "Permutation Bisection";
    worm_off = true;
    worm_on = true;
}

int Permutation_Bisection::attempt(int &id,Parameters &params, Paths &paths, RNG &rng){
    Moves::attempt(id, params, paths, rng);
    if(params.particles == 0) return 0;
    int set = rng.randint(multisteps.size());
    bisection_info.clear();
    bisection_info.push_back(multisteps[set].front());
    bisection_info.push_back(multisteps[set].back());
    if(multisteps[set].back() < multisteps[set].front())
        bisection_info.back() += params.total_slices;
    bisection_info.push_back(rng.randint(params.particles));
    MPI_Bcast(&bisection_info[0], 3, MPI_INT, 0,local_comm);
    if(paths.broken[bisection_info[2]]){
        int final_ptcl = bisection_info[2];
        if(bisection_info[1] > params.total_slices)
            final_ptcl = paths.forward_connects[final_ptcl];
        if((bisection_info[2] == params.worm_head.second && bisection_info[0] < params.worm_head.first) || (bisection_info[2] == params.worm_tail.second && bisection_info[1] > params.worm_tail.first)|| (final_ptcl == params.worm_tail.second && bisection_info[0]+multistep_dist > params.worm_tail.first)){
            return 1;
        }
    }
    keep_going = true;
    int start_slice = bisection_info[0];
    int swap_slice = bisection_info[1]%params.total_slices;
    if(start_slice >= params.my_start && start_slice <= params.my_end){
        if(swap_slice >= params.my_start && swap_slice <= params.my_end){
            std::vector<int> possible_swaps = paths.get_nearest_neighbors(swap_slice%params.slices_per_process, paths.get_coordinate(start_slice%params.slices_per_process, bisection_info[2]));
            std::vector<int> possible_swaps_keys = paths.get_nearest_neighbor_keys(swap_slice%params.slices_per_process, paths.get_coordinate(start_slice%params.slices_per_process, bisection_info[2]));
            if(!possible_swaps.empty())
                for(int poss = possible_swaps.size() - 1; poss >= 0; --poss)
                    if(paths.charge[possible_swaps[poss]] != paths.charge[bisection_info[2]]){
                        possible_swaps.erase(possible_swaps.begin() + poss);
                        possible_swaps_keys.erase(possible_swaps_keys.begin() + poss);
                    }
            if(possible_swaps.size() == 0)
                keep_going = false;
            else{
                std::vector<double> rho0s(possible_swaps.size(),0);
                int sp_fwd = bisection_info[2];
                if(start_slice + params.multistep_dist >= params.total_slices)
                    sp_fwd = paths.forward_connects[bisection_info[2]];
                double r1d = paths.get_kinetic_separation(start_slice%params.slices_per_process, paths.get_coordinate_key(start_slice%params.slices_per_process, bisection_info[2]), paths.get_coordinate_key(swap_slice%params.slices_per_process, sp_fwd));
                for(auto &i : possible_swaps){
                    int ps_back = i;
                    if(start_slice + params.multistep_dist >= params.total_slices)
                        ps_back = paths.backward_connects[i];
                    if(ps_back == -1 || (ps_back == params.worm_head.second && start_slice < params.worm_head.first))
                        rho0s[&i - &possible_swaps[0]] = 0;
                    else{
                        double r1 = paths.get_kinetic_separation(start_slice%params.slices_per_process, paths.get_coordinate_key(start_slice%params.slices_per_process, bisection_info[2]), possible_swaps_keys[&i - &possible_swaps[0]]);
                        double r2 = paths.get_kinetic_separation(start_slice%params.slices_per_process, paths.get_coordinate_key(start_slice%params.slices_per_process, ps_back), paths.get_coordinate_key(swap_slice%params.slices_per_process, sp_fwd));
                        double r2d = paths.get_kinetic_separation(start_slice%params.slices_per_process, paths.get_coordinate_key(start_slice%params.slices_per_process, ps_back), possible_swaps_keys[&i - &possible_swaps[0]]);
                        rho0s[&i - &possible_swaps[0]] = (exp(-1.0/(4.0*params.lambda*params.Mbar*params.tau)*r1)*exp(-1.0/(4.0*params.lambda*params.Mbar*params.tau)*r2))/(exp(-1.0/(4.0*params.lambda*params.Mbar*params.tau)*r1d)*exp(-1.0/(4.0*params.lambda*params.Mbar*params.tau)*r2d));
                    }
                }
                double sum = 0;
                for (auto &r : rho0s)
                    sum += r;
                if(sum < 10E-10)
                    keep_going = false;
                double total = 0;
                double rn = rng.randnormed(sum);
                for(auto &r: rho0s){
                    total += r;
                    if(rn <= total){
                        choice = possible_swaps[&r - &rho0s[0]];
                        break;
                    }
                }
                int choice_back = choice;
                if(start_slice + params.multistep_dist >= params.total_slices)
                    choice_back = paths.backward_connects[choice];
                if(choice_back == bisection_info[2])
                    keep_going = false;
            }
        }
        else{
            MPI_Send(&paths.get_coordinate(start_slice%params.slices_per_process, bisection_info[2])[0], params.dimensions, MPI_DOUBLE,swap_slice/params.slices_per_process, 0, local_comm);
            std::vector<int> possible_swaps;
            std::vector<int> possible_swaps_keys;
            int number_amt = 0;
            MPI_Recv(&number_amt, 1, MPI_INT, swap_slice/params.slices_per_process, 1, local_comm, MPI_STATUS_IGNORE);
            if(number_amt != 0){
                possible_swaps.resize(number_amt);
                possible_swaps_keys.resize(number_amt);
                MPI_Recv(&possible_swaps[0], number_amt, MPI_INT, swap_slice/params.slices_per_process, 2, local_comm, MPI_STATUS_IGNORE);
                MPI_Recv(&possible_swaps_keys[0], number_amt, MPI_INT, swap_slice/params.slices_per_process, 3, local_comm, MPI_STATUS_IGNORE);
                int key_fwd = 0;
                MPI_Recv(&key_fwd, 1, MPI_INT, swap_slice/params.slices_per_process, 4, local_comm, MPI_STATUS_IGNORE);
                std::vector<double> rho0s(possible_swaps.size(),0);
                int sp_fwd = bisection_info[2];
                if(start_slice + params.multistep_dist >= params.total_slices)
                    sp_fwd = paths.forward_connects[bisection_info[2]];
                double r1d = paths.get_kinetic_separation(start_slice%params.slices_per_process, paths.get_coordinate_key(start_slice%params.slices_per_process, bisection_info[2]), key_fwd);
                for(auto &i : possible_swaps){
                    int ps_back = i;
                    if(start_slice + params.multistep_dist >= params.total_slices)
                        ps_back = paths.backward_connects[i];
                    if(ps_back == -1 || (ps_back == params.worm_head.second && start_slice < params.worm_head.first))
                        rho0s[&i - &possible_swaps[0]] = 0;
                    else{
                        double r1 = paths.get_kinetic_separation(start_slice%params.slices_per_process, paths.get_coordinate_key(start_slice%params.slices_per_process, bisection_info[2]), possible_swaps_keys[&i - &possible_swaps[0]]);
                        double r2 = paths.get_kinetic_separation(start_slice%params.slices_per_process, paths.get_coordinate_key(start_slice%params.slices_per_process, ps_back), key_fwd);
                        double r2d = paths.get_kinetic_separation(start_slice%params.slices_per_process, paths.get_coordinate_key(start_slice%params.slices_per_process, ps_back), possible_swaps_keys[&i - &possible_swaps[0]]);
                        rho0s[&i - &possible_swaps[0]] = (exp(-1.0/(4.0*params.lambda*params.Mbar*params.tau)*r1)*exp(-1.0/(4.0*params.lambda*params.Mbar*params.tau)*r2))/(exp(-1.0/(4.0*params.lambda*params.Mbar*params.tau)*r1d)*exp(-1.0/(4.0*params.lambda*params.Mbar*params.tau)*r2d));
                    }
                }
                double sum = 0;
                for (auto &r : rho0s)
                    sum += r;
                if(sum < 10E-10)
                    keep_going = false;
                double total = 0;
                double rn = rng.randnormed(sum);
                for(auto &r: rho0s){
                    total += r;
                    if(rn <= total){
                        choice = possible_swaps[&r - &rho0s[0]];
                        break;
                    }
                }
                int choice_back = choice;
                if(start_slice + params.multistep_dist >= params.total_slices)
                    choice_back = paths.backward_connects[choice];
                if(choice_back == bisection_info[2])
                    keep_going = false;
            }
            else
                keep_going = false;
        }
    }
    else if(swap_slice >= params.my_start && swap_slice <= params.my_end){
        std::vector<double> tail(params.dimensions,0);
        MPI_Recv(&tail[0], params.dimensions, MPI_DOUBLE, start_slice/params.slices_per_process, 0, local_comm, MPI_STATUS_IGNORE);
        int slice = swap_slice%params.slices_per_process;
        std::vector<int> possible_swaps = paths.get_nearest_neighbors(slice, tail);
        std::vector<int> possible_swaps_keys = paths.get_nearest_neighbor_keys(slice, tail);
        if(!possible_swaps.empty())
            for(int poss = possible_swaps.size() - 1; poss >= 0; --poss)
                if(paths.charge[possible_swaps[poss]] != paths.charge[bisection_info[2]]){
                    possible_swaps.erase(possible_swaps.begin() + poss);
                    possible_swaps_keys.erase(possible_swaps_keys.begin() + poss);
                }
        int number_amt = possible_swaps.size();
        MPI_Send(&number_amt, 1, MPI_INT,start_slice/params.slices_per_process, 1, local_comm);
        if(number_amt != 0){
            MPI_Send(&possible_swaps[0], number_amt, MPI_INT, start_slice/params.slices_per_process, 2, local_comm);
            MPI_Send(&possible_swaps_keys[0], number_amt, MPI_INT, start_slice/params.slices_per_process, 3, local_comm);
            int orig_part = bisection_info[2];
            if(bisection_info[1]>=params.total_slices)
                orig_part = paths.forward_connects[bisection_info[2]];
            int key_fwd = paths.get_coordinate_key(slice, orig_part);
            MPI_Send(&key_fwd, 1, MPI_INT,start_slice/params.slices_per_process, 4, local_comm);
        }
    }
    MPI_Bcast(&keep_going, 1, MPI_INT, start_slice/params.slices_per_process, local_comm);
    if(!keep_going)
        return 1;
    MPI_Bcast(&choice, 1, MPI_INT, start_slice/params.slices_per_process, local_comm);
    new_coordinates.resize(params.slices_per_process*2);
    new_distances.reserve(params.slices_per_process*params.particles*2);
    new_coordinates_ahead.resize(params.slices_per_process*2);
    keys_ahead.resize(params.slices_per_process*2);
    int level = int(round(log2(bisection_info[1]-bisection_info[0])));
    int counter = 0;
    ptcl_slice_1.resize(params.slices_per_process,bisection_info[2]);
    ptcl_slice_2.resize(params.slices_per_process,choice);
    if(bisection_info[1]>= params.total_slices)
        for(int i = bisection_info[0]; i <=bisection_info[1]; ++i){
            if(i%params.total_slices >= params.my_start && i%params.total_slices <= params.my_end && i >= params.total_slices)
                ptcl_slice_1[i%params.slices_per_process] = paths.forward_connects[bisection_info[2]];
            if(i%params.total_slices >= params.my_start && i%params.total_slices <= params.my_end && i < params.total_slices)
                ptcl_slice_2[i%params.slices_per_process] = paths.backward_connects[choice];
        }
    
    while(level){
        int num_comps = pow(2,counter);
        counter++;
        std::vector<int> slices_to_compute;
        for(int i = 0; i < num_comps; ++i)
            slices_to_compute.push_back(int((1/pow(2.,counter) + i/pow(2.,counter-1))*(bisection_info[1]-bisection_info[0])+bisection_info[0])%params.total_slices);
        for(auto &slice : slices_to_compute){
            int first = int(slice - pow(2,level-1)+params.total_slices)%params.total_slices;
            int second = int(slice + pow(2, level-1))%params.total_slices;
            if(slice >= params.my_start && slice <= params.my_end){
                std::vector<double> start_1(params.dimensions);
                std::vector<double> end_1(params.dimensions);
                std::vector<double> start_2(params.dimensions);
                std::vector<double> end_2(params.dimensions);
                if(first < params.my_start || first > params.my_end){
                    MPI_Recv(&start_1[0], params.dimensions,MPI_DOUBLE,first/params.slices_per_process,2*counter-1,local_comm,MPI_STATUS_IGNORE);
                    MPI_Recv(&start_2[0], params.dimensions,MPI_DOUBLE,first/params.slices_per_process,2*counter-1,local_comm,MPI_STATUS_IGNORE);
                }
                else{
                    if(counter == 1){
                        start_1 = paths.get_coordinate(first%params.slices_per_process, ptcl_slice_1[first%params.slices_per_process]);
                        start_2 = paths.get_coordinate(first%params.slices_per_process, ptcl_slice_2[first%params.slices_per_process]);
                        new_coordinates[first%params.slices_per_process] = start_1;
                        new_coordinates[first%params.slices_per_process+params.slices_per_process] = start_2;
                    }
                    else{
                        start_1 = new_coordinates[first%params.slices_per_process];
                        start_2 = new_coordinates[first%params.slices_per_process+params.slices_per_process];
                    }
                }
                if(second < params.my_start || second > params.my_end){
                    MPI_Recv(&end_1[0], params.dimensions,MPI_DOUBLE,second/params.slices_per_process,2*counter,local_comm,MPI_STATUS_IGNORE);
                    MPI_Recv(&end_2[0], params.dimensions,MPI_DOUBLE,second/params.slices_per_process,2*counter,local_comm,MPI_STATUS_IGNORE);
                }
                else{
                    if(counter == 1){
                        end_1 = paths.get_coordinate(second%params.slices_per_process, ptcl_slice_2[second%params.slices_per_process]);
                        end_2 = paths.get_coordinate(second%params.slices_per_process, ptcl_slice_1[second%params.slices_per_process]);
                        new_coordinates[second%params.slices_per_process] = end_1;
                        new_coordinates[second%params.slices_per_process+params.slices_per_process] = end_2;
                    }
                    else{
                        end_1 = new_coordinates[second%params.slices_per_process];
                        end_2 = new_coordinates[second%params.slices_per_process+params.slices_per_process];
                    }
                }
                double tau = params.tau*pow(2, level-1);
                double width = sqrt(params.lambda*tau);
                average_loc(start_1, end_1, new_coordinates[slice%params.slices_per_process], params.box_size);
                for(auto &dim : new_coordinates[slice%params.slices_per_process]){
                    dim += rng.randgaussian(width);
                }
                average_loc(start_2, end_2, new_coordinates[slice%params.slices_per_process+params.slices_per_process], params.box_size);
                for(auto &dim : new_coordinates[slice%params.slices_per_process+params.slices_per_process]){
                    dim += rng.randgaussian(width);
                }
            }
            else{
                if(first >= params.my_start && first <= params.my_end){
                    if(counter == 1){
                        new_coordinates[first%params.slices_per_process] = paths.get_coordinate(first%params.slices_per_process, ptcl_slice_1[first%params.slices_per_process]);
                        new_coordinates[first%params.slices_per_process+params.slices_per_process] = paths.get_coordinate(first%params.slices_per_process, ptcl_slice_2[first%params.slices_per_process]);
                    }
                    MPI_Send(&new_coordinates[first%params.slices_per_process][0],params.dimensions,MPI_DOUBLE,slice/params.slices_per_process,2*counter-1,local_comm);
                    MPI_Send(&new_coordinates[first%params.slices_per_process+params.slices_per_process][0],params.dimensions,MPI_DOUBLE,slice/params.slices_per_process,2*counter-1,local_comm);
                }
                if(second >= params.my_start && second <= params.my_end){
                    if(counter == 1){
                        new_coordinates[second%params.slices_per_process] = paths.get_coordinate(second%params.slices_per_process, ptcl_slice_2[second%params.slices_per_process]);
                        new_coordinates[second%params.slices_per_process+params.slices_per_process] = paths.get_coordinate(second%params.slices_per_process, ptcl_slice_1[second%params.slices_per_process]);
                    }
                    MPI_Send(&new_coordinates[second%params.slices_per_process][0],params.dimensions,MPI_DOUBLE,slice/params.slices_per_process,2*counter,local_comm);
                    MPI_Send(&new_coordinates[second%params.slices_per_process+params.slices_per_process][0],params.dimensions,MPI_DOUBLE,slice/params.slices_per_process,2*counter,local_comm);
                }
            }
        }
        --level;
    }
    check(id, params, paths, rng);
    MPI_Bcast(&ac_re, 1, MPI_INT, 0, local_comm);
    ++num_attempts;
    if(ac_re){
        ++num_accepts;
        int p1 = bisection_info[2];
        int p2 = choice;
        if(bisection_info[1] >= params.total_slices)
            p2 = paths.backward_connects[choice];
        paths.swap_worldlines(params, bisection_info[0], p1, p2);
        int sp1 = p1;
        int sp2 = p2;
        for(int i = bisection_info[0]; i <= bisection_info[1]; ++i){
            if(i == params.total_slices){
                p1 = paths.forward_connects[p1];
                p2 = paths.forward_connects[p2];
            }
            if(i%params.total_slices >= params.my_start && i%params.total_slices <= params.my_end){
                int slice = (i%params.total_slices)%params.slices_per_process;
                put_in_box(new_coordinates[slice], params.box_size);
                paths.set_coordinate(slice, p1, new_coordinates[slice], false);
                put_in_box(new_coordinates[slice+params.slices_per_process], params.box_size);
                paths.set_coordinate(slice, p2, new_coordinates[slice+params.slices_per_process], false);
            }
        }
        paths.update_separations(new_distances);
        if(!params.gce){
            p1 = sp1;
            p2 = sp2;
            for(int i = bisection_info[0]; i <= bisection_info[1]; ++i){
                if(i == params.total_slices){
                    p1 = paths.forward_connects[p1];
                    p2 = paths.forward_connects[p2];
                }
                if(i%params.total_slices >= params.my_start && i%params.total_slices <= params.my_end){
                    int slice = (i%params.total_slices)%params.slices_per_process;
                    paths.calculate_kinetic_separations_start(slice, paths.get_coordinate_key(slice, p1));
                    paths.calculate_kinetic_separations_start(slice, paths.get_coordinate_key(slice, p2));
                }
            }
            p1 = sp1;
            p2 = sp2;
            for(int i = bisection_info[0]; i <= bisection_info[1]; ++i){
                if(i == params.total_slices){
                    p1 = paths.forward_connects[p1];
                    p2 = paths.forward_connects[p2];
                }
                int slice_back = positive_modulo(i-params.multistep_dist, params.total_slices);
                int tslice = i%params.total_slices;
                if(tslice >= params.my_start && tslice <= params.my_end){
                    int slice = tslice%params.slices_per_process;
                    int send_rank = slice_back/params.slices_per_process;
                    int send_tag = slice_back;
                    new_coordinates[slice].resize(params.dimensions);
                    int key_1 = paths.get_coordinate_key(slice, p1);
                    int key_2 = paths.get_coordinate_key(slice, p2);
                    if(send_rank != id){
                        MPI_Send(&new_coordinates[slice][0], params.dimensions, MPI_DOUBLE, send_rank, send_tag, local_comm);
                        MPI_Send(&key_1, 1, MPI_INT, send_rank, send_tag+params.slices_per_process, local_comm);
                        MPI_Send(&new_coordinates[slice+params.slices_per_process][0], params.dimensions, MPI_DOUBLE, send_rank, send_tag+2*params.slices_per_process, local_comm);
                        MPI_Send(&key_2, 1, MPI_INT, send_rank, send_tag+3*params.slices_per_process, local_comm);
                    }
                    else{
                        keys_ahead[slice_back%params.slices_per_process] = key_1;
                        keys_ahead[slice_back%params.slices_per_process+params.slices_per_process] = key_2;
                        new_coordinates_ahead[slice_back%params.slices_per_process] = new_coordinates[slice];
                        new_coordinates_ahead[slice_back%params.slices_per_process+params.slices_per_process] = new_coordinates[slice+params.slices_per_process];
                        paths.set_kinetic_end(slice_back%params.slices_per_process, keys_ahead[slice_back%params.slices_per_process], new_coordinates_ahead[slice_back%params.slices_per_process], true);
                        paths.set_kinetic_end(slice_back%params.slices_per_process, keys_ahead[slice_back%params.slices_per_process+params.slices_per_process], new_coordinates_ahead[slice_back%params.slices_per_process+params.slices_per_process], true);
                    }
                }
                else if(slice_back >= params.my_start && slice_back <= params.my_end){
                    int slice = slice_back%params.slices_per_process;
                    int recv_rank = tslice/params.slices_per_process;
                    int recv_tag = slice_back;
                    new_coordinates_ahead[slice].resize(params.dimensions);
                    new_coordinates_ahead[slice+params.slices_per_process].resize(params.dimensions);
                    MPI_Recv(&new_coordinates_ahead[slice][0], params.dimensions, MPI_DOUBLE, recv_rank, recv_tag, local_comm, MPI_STATUS_IGNORE);
                    MPI_Recv(&keys_ahead[slice], 1, MPI_INT, recv_rank, recv_tag+params.slices_per_process, local_comm, MPI_STATUS_IGNORE);
                    MPI_Recv(&new_coordinates_ahead[slice+params.slices_per_process][0], params.dimensions, MPI_DOUBLE, recv_rank, recv_tag+2*params.slices_per_process, local_comm, MPI_STATUS_IGNORE);
                    MPI_Recv(&keys_ahead[slice+params.slices_per_process], 1, MPI_INT, recv_rank, recv_tag+3*params.slices_per_process, local_comm, MPI_STATUS_IGNORE);
                    paths.set_kinetic_end(slice, keys_ahead[slice], new_coordinates_ahead[slice], true);
                    paths.set_kinetic_end(slice, keys_ahead[slice+params.slices_per_process], new_coordinates_ahead[slice+params.slices_per_process], true);
                }
            }
        }
    }
    bisection_info.clear();
    ptcl_slice_1.clear();
    ptcl_slice_2.clear();
    new_coordinates.clear();
    new_coordinates_ahead.clear();
    keys_ahead.clear();
    new_distances.clear();
    return 0;
}

void Permutation_Bisection::check(int &id, Parameters &params, Paths &paths, RNG &rng){
    std::vector<double> old_action_p(params.slices_per_process,0);
    std::vector<double> new_action_p(params.slices_per_process,0);
    std::vector<double> dist(params.dimensions);
    std::vector<int> kVec(params.dimensions);
    if(params.potential == 2){
        params.alpha = sqrt(M_PI)*pow(params.particles/params.volume2,1/6.);
        params.coulcut = sqrt(params.p)/params.alpha;
        params.coulcut2 = pow(params.coulcut,2);
        params.kcut = 2.*params.p/params.coulcut;
        params.kcut2 = pow(params.kcut,2);
        params.nmax = floor(params.coulcut/params.box_size);
        params.kmax = ceil(params.kcut/(2.*M_PI/params.box_size));
    }
    for(int i = bisection_info[0]; i <= bisection_info[1]; ++i){
        int slice =i%params.total_slices;
        if(slice >= params.my_start && slice <= params.my_end){
            double old_action_local = 0;
            double new_action_local = 0;
            switch(params.potential){
                case 0:
                    old_action_p[slice%params.slices_per_process] += potential_value_harmonic(inner_product(paths.get_coordinate(slice%params.slices_per_process,ptcl_slice_1[slice%params.slices_per_process]).begin(),paths.get_coordinate(slice%params.slices_per_process,ptcl_slice_1[slice%params.slices_per_process]).end(),paths.get_coordinate(slice%params.slices_per_process,ptcl_slice_1[slice%params.slices_per_process]).begin(),0.0));
                    old_action_p[slice%params.slices_per_process] += potential_value_harmonic(inner_product(paths.get_coordinate(slice%params.slices_per_process,ptcl_slice_2[slice%params.slices_per_process]).begin(),paths.get_coordinate(slice%params.slices_per_process,ptcl_slice_2[slice%params.slices_per_process]).end(),paths.get_coordinate(slice%params.slices_per_process,ptcl_slice_2[slice%params.slices_per_process]).begin(),0.0));
                    new_action_p[slice%params.slices_per_process] += potential_value_harmonic(inner_product(new_coordinates[slice%params.slices_per_process].begin(),new_coordinates[slice%params.slices_per_process].end(),new_coordinates[slice%params.slices_per_process].begin(),0.0));
                    new_action_p[slice%params.slices_per_process] += potential_value_harmonic(inner_product(new_coordinates[slice%params.slices_per_process+params.slices_per_process].begin(),new_coordinates[slice%params.slices_per_process+params.slices_per_process].end(),new_coordinates[slice%params.slices_per_process+params.slices_per_process].begin(),0.0));
                    break;
                case 1:
                {
#pragma omp parallel for firstprivate(dist) reduction (+:new_action_local, old_action_local)
                    for(int ptcl = 0; ptcl < params.particles; ++ptcl){
                        if(!(ptcl == params.worm_head.second && slice < params.worm_head.first) && !(ptcl == params.worm_tail.second && slice > params.worm_tail.first)){
                            if(ptcl != ptcl_slice_1[slice%params.slices_per_process] && ptcl != ptcl_slice_2[slice%params.slices_per_process]){
                                old_action_local += potential_value_aziz(paths.get_separation(slice%params.slices_per_process, ptcl, ptcl_slice_1[slice%params.slices_per_process]));
                                old_action_local += potential_value_aziz(paths.get_separation(slice%params.slices_per_process, ptcl, ptcl_slice_2[slice%params.slices_per_process]));
                                distance(paths.get_coordinate(slice%params.slices_per_process,ptcl), new_coordinates[slice%params.slices_per_process], dist, params.box_size);
                                double r = sqrt(std::inner_product(dist.begin(), dist.end(), dist.begin(), 0.0));
                                new_action_local += potential_value_aziz(r);
#pragma omp critical
                                new_distances.push_back(std::tuple<std::pair<int,int>,std::vector<double>,double>(std::pair<int,int>(paths.get_coordinate_key(slice%params.slices_per_process, ptcl),paths.get_coordinate_key(slice%params.slices_per_process, ptcl_slice_2[slice%params.slices_per_process])), dist, r));
                                distance(paths.get_coordinate(slice%params.slices_per_process,ptcl), new_coordinates[slice%params.slices_per_process+params.slices_per_process], dist, params.box_size);
                                r = sqrt(std::inner_product(dist.begin(), dist.end(), dist.begin(), 0.0));
                                new_action_local += potential_value_aziz(r);
#pragma omp critical
                                new_distances.push_back(std::tuple<std::pair<int,int>,std::vector<double>,double>(std::pair<int,int>(paths.get_coordinate_key(slice%params.slices_per_process, ptcl),paths.get_coordinate_key(slice%params.slices_per_process, ptcl_slice_1[slice%params.slices_per_process])), dist, r));
                            }
                        }
                    }
                    old_action_local += 2*potential_value_aziz(paths.get_separation(slice%params.slices_per_process, ptcl_slice_2[slice%params.slices_per_process], ptcl_slice_1[slice%params.slices_per_process]));
                    distance(new_coordinates[slice%params.slices_per_process+params.slices_per_process], new_coordinates[slice%params.slices_per_process], dist, params.box_size);
                    double r = sqrt(std::inner_product(dist.begin(), dist.end(), dist.begin(), 0.0));
                    new_action_local += 2*potential_value_aziz(r);
                    new_distances.push_back(std::tuple<std::pair<int,int>,std::vector<double>,double>(std::pair<int,int>(paths.get_coordinate_key(slice%params.slices_per_process, ptcl_slice_1[slice%params.slices_per_process]),paths.get_coordinate_key(slice%params.slices_per_process, ptcl_slice_2[slice%params.slices_per_process])), dist, r));
                    old_action_p[slice%params.slices_per_process] += old_action_local;
                    new_action_p[slice%params.slices_per_process] += new_action_local;
                    break;
                }
                case 2:
                {
#pragma omp parallel for firstprivate(dist) reduction (+:new_action_local, old_action_local)
                    for(int ptcl = 0; ptcl < params.particles; ++ptcl){
                        if(!(ptcl == params.worm_head.second && slice < params.worm_head.first) && !(ptcl == params.worm_tail.second && slice > params.worm_tail.first)){
                            if(ptcl != ptcl_slice_1[slice%params.slices_per_process] && ptcl != ptcl_slice_2[slice%params.slices_per_process]){
                                old_action_local += potential_value_coulomb(paths.get_separation_vector(slice%params.slices_per_process, ptcl, ptcl_slice_1[slice%params.slices_per_process]), kVec,  paths.charge[ptcl],paths.charge[ptcl_slice_1[slice%params.slices_per_process]], params);
                                old_action_local += potential_value_coulomb(paths.get_separation_vector(slice%params.slices_per_process, ptcl, ptcl_slice_2[slice%params.slices_per_process]), kVec,  paths.charge[ptcl],paths.charge[ptcl_slice_2[slice%params.slices_per_process]], params);
                                distance(paths.get_coordinate(slice%params.slices_per_process,ptcl), new_coordinates[slice%params.slices_per_process], dist, params.box_size);
                                double r = sqrt(std::inner_product(dist.begin(), dist.end(), dist.begin(), 0.0));
                                new_action_local += potential_value_coulomb(dist, kVec, paths.charge[ptcl], paths.charge[ptcl_slice_1[slice%params.slices_per_process]],  params);
#pragma omp critical
                                new_distances.push_back(std::tuple<std::pair<int,int>,std::vector<double>,double>(std::pair<int,int>(paths.get_coordinate_key(slice%params.slices_per_process, ptcl),paths.get_coordinate_key(slice%params.slices_per_process, ptcl_slice_2[slice%params.slices_per_process])), dist, r));
                                distance(paths.get_coordinate(slice%params.slices_per_process,ptcl), new_coordinates[slice%params.slices_per_process+params.slices_per_process], dist, params.box_size);
                                r = sqrt(std::inner_product(dist.begin(), dist.end(), dist.begin(), 0.0));
                                new_action_local += potential_value_coulomb(dist, kVec, paths.charge[ptcl], paths.charge[ptcl_slice_2[slice%params.slices_per_process]],  params);
#pragma omp critical
                                new_distances.push_back(std::tuple<std::pair<int,int>,std::vector<double>,double>(std::pair<int,int>(paths.get_coordinate_key(slice%params.slices_per_process, ptcl),paths.get_coordinate_key(slice%params.slices_per_process, ptcl_slice_1[slice%params.slices_per_process])), dist, r));

                            }
                        }
                    }
                    old_action_local += 2*potential_value_coulomb(paths.get_separation_vector(slice%params.slices_per_process, ptcl_slice_2[slice%params.slices_per_process], ptcl_slice_1[slice%params.slices_per_process]), kVec,  paths.charge[ptcl_slice_2[slice%params.slices_per_process]],paths.charge[ptcl_slice_1[slice%params.slices_per_process]], params);
                    distance(new_coordinates[slice%params.slices_per_process+params.slices_per_process], new_coordinates[slice%params.slices_per_process], dist, params.box_size);
                    double r = sqrt(std::inner_product(dist.begin(), dist.end(), dist.begin(), 0.0));
                    new_action_local += 2*potential_value_coulomb(dist, kVec,paths.charge[ptcl_slice_2[slice%params.slices_per_process]],paths.charge[ptcl_slice_1[slice%params.slices_per_process]], params);
                    new_distances.push_back(std::tuple<std::pair<int,int>,std::vector<double>,double>(std::pair<int,int>(paths.get_coordinate_key(slice%params.slices_per_process, ptcl_slice_1[slice%params.slices_per_process]),paths.get_coordinate_key(slice%params.slices_per_process, ptcl_slice_2[slice%params.slices_per_process])), dist, r));
                    old_action_p[slice%params.slices_per_process] = old_action_local;
                    new_action_p[slice%params.slices_per_process] = new_action_local;
                    break;
                }
            }
        }
    }
    if(id != 0){
        MPI_Send(&old_action_p[0], params.slices_per_process, MPI_DOUBLE, 0, 0, local_comm);
        MPI_Send(&new_action_p[0], params.slices_per_process, MPI_DOUBLE, 0, 1, local_comm);
    }
    else{
        ac_re = false;
        std::vector<double> old_action;
        std::vector<double> new_action;
        old_action.reserve(params.total_slices);
        new_action.reserve(params.total_slices);
        old_action.insert(old_action.end(),old_action_p.begin(),old_action_p.end());
        new_action.insert(new_action.end(),new_action_p.begin(),new_action_p.end());
        for(int i = 1; i < params.num_workers; ++i){
            MPI_Recv(&old_action_p[0], params.slices_per_process, MPI_DOUBLE, i, 0, local_comm, MPI_STATUS_IGNORE);
            MPI_Recv(&new_action_p[0], params.slices_per_process, MPI_DOUBLE, i, 1, local_comm, MPI_STATUS_IGNORE);
            old_action.insert(old_action.end(),old_action_p.begin(),old_action_p.end());
            new_action.insert(new_action.end(),new_action_p.begin(),new_action_p.end());
        }
        double oa = 0;
        double na = 0;
        for(int j = bisection_info[0]; j < bisection_info[1]; ++j){
            int slice = j%params.total_slices;
            int slicep1 = (j+1)%params.total_slices;
            oa += params.tau/2.*(old_action[slice] + old_action[slicep1]);
            na += params.tau/2.*(new_action[slice] + new_action[slicep1]);
        }
        if(rng.randnormed(1) < exp(-(na - oa)))
            ac_re = true;
    }
}

int Open::attempt(int &id, Parameters &params, Paths &paths, RNG &rng){
    if(params.worm_on) return 0;
    if(params.particles == 0) return 0;
    Moves::attempt(id, params, paths, rng);    open_info.push_back(rng.randint(params.particles));
    open_info.push_back(rng.randint(params.total_slices));
    open_info.push_back(rng.randint(params.Mbar)+1);
    MPI_Bcast(&open_info[0], 3, MPI_INT, 0, local_comm);
    check(id, params, paths, rng);
    MPI_Bcast(&ac_re, 1, MPI_INT,0, local_comm);
    ++num_attempts;
    if(open_info[1]-open_info[2] < 0){
        open_info[0] = paths.backward_connects[open_info[0]];
    }
    if(ac_re){
        ++num_accepts;
        paths.open_path(params,open_info[0], positive_modulo(open_info[1]-open_info[2],params.total_slices),open_info[2]);
    }
    open_info.clear();
    return 0;
}

void Open::check(int &id, Parameters &params, Paths &paths, RNG &rng){
    std::vector<double> first_part(params.dimensions,0);
    std::vector<double> second_part(params.dimensions,0);
    std::vector<int> kVec(params.dimensions);
    int ss = open_info[1]-open_info[2];
    int es = open_info[1];
    int myp = open_info[0];
    std::vector<double> old_action_p(params.slices_per_process,0);
    if(params.potential == 2){
        params.alpha = sqrt(M_PI)*pow(params.particles/params.volume2,1/6.);
        params.coulcut = sqrt(params.p)/params.alpha;
        params.coulcut2 = pow(params.coulcut,2);
        params.kcut = 2.*params.p/params.coulcut;
        params.kcut2 = pow(params.kcut,2);
        params.nmax = floor(params.coulcut/params.box_size);
        params.kmax = ceil(params.kcut/(2.*M_PI/params.box_size));
    }
    for(int i = open_info[1]-open_info[2]+1; i < open_info[1]; ++i){
        int slice = (i+params.total_slices)%params.total_slices;
        int part = open_info[0];
        if(i < 0) part = paths.backward_connects[open_info[0]];
        if(slice >= params.my_start && slice <= params.my_end){
            double old_action_local = 0;
            switch(params.potential){
                case 0:
                    old_action_p[slice%params.slices_per_process] = potential_value_harmonic(inner_product(paths.get_coordinate(slice%params.slices_per_process,part).begin(),paths.get_coordinate(slice%params.slices_per_process,part).end(),paths.get_coordinate(slice%params.slices_per_process,part).begin(),0.0));
                    break;
                case 1:
#pragma omp parallel for reduction (+:old_action_local)
                    for(int ptcl = 0; ptcl < params.particles; ++ptcl){
                        if(!(ptcl == params.worm_head.second && slice < params.worm_head.first) && !(ptcl == params.worm_tail.second && slice > params.worm_tail.first)){
                            if(ptcl != part){
                                 old_action_local += potential_value_aziz(paths.get_separation(slice%params.slices_per_process, ptcl, part));
                            }
                        }
                    }
                    old_action_p[slice%params.slices_per_process] += old_action_local;
                    break;
                case 2:
#pragma omp parallel for reduction (+:old_action_local)
                    for(int ptcl = 0; ptcl < params.particles; ++ptcl){
                        if(!(ptcl == params.worm_head.second && slice < params.worm_head.first) && !(ptcl == params.worm_tail.second && slice > params.worm_tail.first)){
                            old_action_local += potential_value_coulomb(paths.get_separation_vector(slice%params.slices_per_process,  ptcl, part), kVec, paths.charge[ptcl], paths.charge[part], params);
                        }
                    }
                    old_action_p[slice%params.slices_per_process] += old_action_local;
                    break;
            }
        }
    }
    if(id != 0){
        if((ss+params.total_slices)%params.total_slices >= params.my_start && (ss+params.total_slices)%params.total_slices <= params.my_end){
            if(ss < 0){
                ss = ss+params.total_slices;
                myp = paths.backward_connects[open_info[0]];
            }
            first_part = paths.get_coordinate(ss%params.slices_per_process,myp);
            MPI_Send(&first_part[0],params.dimensions,MPI_DOUBLE,0,1,local_comm);
        }
        if(es >= params.my_start && es <= params.my_end){
            second_part = paths.get_coordinate(es%params.slices_per_process,myp);
            MPI_Send(&second_part[0],params.dimensions,MPI_DOUBLE,0,2,local_comm);
        }
        MPI_Send(&old_action_p[0], params.slices_per_process, MPI_DOUBLE, 0,0,local_comm);
    }
    else{
        MPI_Request req[2];
        if((ss+params.total_slices)%params.total_slices >= params.my_start && (ss+params.total_slices)%params.total_slices <= params.my_end){
            if(ss < 0){
                ss = ss+params.total_slices;
                myp = paths.backward_connects[open_info[0]];
            }
            first_part = paths.get_coordinate(ss%params.slices_per_process,myp);
            req[0] = MPI_REQUEST_NULL;
        }
        else
            MPI_Irecv(&first_part[0],params.dimensions,MPI_DOUBLE,((open_info[1]-open_info[2]+params.total_slices)%params.total_slices)/params.slices_per_process,1,local_comm, &req[0]);
        if(es >= params.my_start && es <= params.my_end){
            second_part = paths.get_coordinate(es%params.slices_per_process,myp);
            req[1] = MPI_REQUEST_NULL;
        }
        else
            MPI_Irecv(&second_part[0],params.dimensions,MPI_DOUBLE,open_info[1]/params.slices_per_process,2,local_comm, &req[1]);
        MPI_Waitall(2,req, MPI_STATUS_IGNORE);
        ac_re = false;
        std::vector<double> old_action;
        old_action.reserve(params.total_slices);
        old_action.insert(old_action.end(),old_action_p.begin(),old_action_p.end());
        for(int j = 1; j < params.num_workers; ++j){
            MPI_Recv(&old_action_p[0], params.slices_per_process, MPI_DOUBLE,j,0,local_comm,MPI_STATUS_IGNORE);
            old_action.insert(old_action.end(),old_action_p.begin(),old_action_p.end());
        }
        std::vector<double> dist;
        distance(first_part, second_part, dist, params.box_size);
        double r2 = inner_product(dist.begin(),dist.end(), dist.begin(),0.0);
        double rho0 = pow((4* M_PI*params.lambda*open_info[2]*params.tau),-params.dimensions/2.)*exp(-1.0/(4.0*params.lambda*open_info[2]*params.tau)*r2);
        double oa = 0;
        for(int i = open_info[1]-open_info[2]; i < open_info[1]; ++i)
            oa += params.tau/2*(old_action[(i+params.total_slices)%params.total_slices]+old_action[(i+1+params.total_slices)%params.total_slices]);
        if(rng.randnormed(1) < (params.C*params.Mbar*params.particles*params.total_slices*exp(oa)*exp(-params.mu*open_info[2]*params.tau))/rho0)
            ac_re = true;
    }
}

int Close::attempt(int &id, Parameters &params, Paths &paths, RNG &rng){
    if(!params.worm_on) return 0;
    int distance = (params.worm_head.first - params.worm_tail.first + params.total_slices)%params.total_slices;
    if(distance > params.Mbar + 1 || distance == 0) return 1;
    Moves::attempt(id, params, paths, rng);
    new_coordinates.resize(params.slices_per_process);
    new_distances.resize(params.slices_per_process);
    std::vector<double> end(params.dimensions,0);
    std::vector<double> start(params.dimensions,0);
    if(params.worm_head.first >= params.my_start && params.worm_head.first <= params.my_end){
        end = paths.get_coordinate(params.worm_head.first%params.slices_per_process, params.worm_head.second);
        new_coordinates[params.worm_head.first%params.slices_per_process] = end;
    }
    MPI_Bcast(&end[0],params.dimensions,MPI_DOUBLE,params.worm_head.first/params.slices_per_process,local_comm);
    int start_slice = (params.worm_tail.first+1)%params.total_slices;
    int start_col = params.worm_tail.second;
    for(int i = 0; i < distance - 1; ++i){
        int slice = (start_slice+i)%params.total_slices;
        int slicem1 = (start_slice+i-1+params.total_slices)%params.total_slices;
        if(slice >= params.my_start && slice <= params.my_end){
            if(slicem1 >= params.my_start && slicem1 <= params.my_end){
                if(i == 0){
                    new_coordinates[slicem1%params.slices_per_process] = paths.get_coordinate(slicem1%params.slices_per_process,start_col);
                }
                start = new_coordinates[slicem1%params.slices_per_process];
            }
            else
                MPI_Recv(&start[0], params.dimensions, MPI_DOUBLE, slicem1/params.slices_per_process, i, local_comm, MPI_STATUS_IGNORE);
            average_loc_weighted(start, end, new_coordinates[slice%params.slices_per_process], params.box_size, 1, distance-1-i);
            double tau_w = params.tau*distance-1-i;
            double width = sqrt(2.*params.lambda/(1./params.tau+1./tau_w));
            for(auto &j : new_coordinates[slice%params.slices_per_process])
                j += rng.randgaussian(width);
            put_in_box(new_coordinates[slice%params.slices_per_process], params.box_size);
        }
        else{
            if(slicem1 >= params.my_start && slicem1 <= params.my_end){
                if(i == 0){
                    new_coordinates[slicem1%params.slices_per_process] = paths.get_coordinate(slicem1%params.slices_per_process,start_col);
                }
                MPI_Send(&new_coordinates[slicem1%params.slices_per_process][0], params.dimensions,MPI_DOUBLE,slice/params.slices_per_process, i, local_comm);
            }
        }
    }
    if(distance == 1 && params.worm_tail.first >= params.my_start && params.worm_tail.first <= params.my_end)
        new_coordinates[params.worm_tail.first%params.slices_per_process] = paths.get_coordinate(params.worm_tail.first%params.slices_per_process,start_col);
    check(id, params, paths, rng);
    MPI_Bcast(&ac_re, 1, MPI_INT,0, local_comm);
    ++num_attempts;
    if(ac_re){
        ++num_accepts;
        paths.close_worm(params, new_coordinates, new_distances);
    }
    new_coordinates.clear();
    new_distances.clear();
    return 0;
}

void Close::check(int &id, Parameters &params, Paths &paths, RNG &rng){
    int disty = (params.worm_head.first - params.worm_tail.first + params.total_slices)%params.total_slices;
    std::vector<double> first_part(params.dimensions,0);
    std::vector<double> second_part(params.dimensions,0);
    std::vector<double> new_action_p(params.slices_per_process,0);
    std::vector<double> dist(params.dimensions,0);
    std::vector<int> kVec(params.dimensions);
    if(params.potential == 2){
        params.alpha = sqrt(M_PI)*pow(params.particles/params.volume2,1/6.);
        params.coulcut = sqrt(params.p)/params.alpha;
        params.coulcut2 = pow(params.coulcut,2);
        params.kcut = 2.*params.p/params.coulcut;
        params.kcut2 = pow(params.kcut,2);
        params.nmax = floor(params.coulcut/params.box_size);
        params.kmax = ceil(params.kcut/(2.*M_PI/params.box_size));
    }
    int start_slice = (params.worm_tail.first)%params.total_slices;
    for(int i = 0; i <= disty; ++i){
        int slice = (start_slice+i)%params.total_slices;
        if(slice >= params.my_start && slice <= params.my_end){
            double new_action_local = 0.0;
            new_distances[slice%params.slices_per_process].reserve(params.particles);
            switch(params.potential){
                case 0:
                    new_action_p[slice%params.slices_per_process] = potential_value_harmonic(inner_product(new_coordinates[slice%params.slices_per_process].begin(), new_coordinates[slice%params.slices_per_process].end(), new_coordinates[slice%params.slices_per_process].begin(), 0.0));
                    break;
                case 1:
#pragma omp parallel for firstprivate(dist) reduction (+:new_action_local)
                    for(int ptcl = 0; ptcl < params.particles; ++ptcl){
                        if(!(ptcl == params.worm_head.second && slice <= params.worm_head.first) && !(ptcl == params.worm_tail.second && slice >= params.worm_tail.first)){
                            distance(paths.get_coordinate(slice%params.slices_per_process,ptcl), new_coordinates[slice%params.slices_per_process], dist, params.box_size);
                            double r = sqrt(std::inner_product(dist.begin(), dist.end(), dist.begin(), 0.0));
                            new_action_local += potential_value_aziz(r);
#pragma omp critical
                            new_distances[slice%params.slices_per_process].push_back(std::tuple<int,std::vector<double>,double>(paths.get_coordinate_key(slice%params.slices_per_process, ptcl), dist, r));
                        }
                    }
                    new_action_p[slice%params.slices_per_process] += new_action_local;
                    break;
                case 2:
#pragma omp parallel for firstprivate(dist) reduction (+:new_action_local)
                    for(int ptcl = 0; ptcl < params.particles; ++ptcl){
                        if(!(ptcl == params.worm_head.second && slice <= params.worm_head.first) && !(ptcl == params.worm_tail.second && slice >= params.worm_tail.first)){
                            distance(paths.get_coordinate(slice%params.slices_per_process,ptcl), new_coordinates[slice%params.slices_per_process], dist, params.box_size);
                            double r = sqrt(std::inner_product(dist.begin(), dist.end(), dist.begin(), 0.0));
                            new_action_local += potential_value_coulomb(dist, kVec, paths.charge[ptcl],paths.charge[params.worm_tail.second], params);
#pragma omp critical
                            new_distances[slice%params.slices_per_process].push_back(std::tuple<int,std::vector<double>,double>(paths.get_coordinate_key(slice%params.slices_per_process, ptcl), dist, r));
                        }
                    }
                    new_action_p[slice%params.slices_per_process] += new_action_local;
                    dist = std::vector<double>(params.dimensions, 0);
                    new_action_p[slice%params.slices_per_process] += potential_value_coulomb(dist, kVec, paths.charge[params.worm_tail.second],paths.charge[params.worm_tail.second], params);
                    break;

            }
        }
    }
    if(id != 0){
        if(params.worm_tail.first >= params.my_start && params.worm_tail.first <= params.my_end)
            MPI_Send(&new_coordinates[params.worm_tail.first%params.slices_per_process][0],params.dimensions,MPI_DOUBLE,0,1,local_comm);
        if(params.worm_head.first >= params.my_start && params.worm_head.first <= params.my_end)
            MPI_Send(&new_coordinates[params.worm_head.first%params.slices_per_process][0],params.dimensions,MPI_DOUBLE,0,2,local_comm);
        MPI_Send(&new_action_p[0], params.slices_per_process, MPI_DOUBLE, 0,0,local_comm);
    }
    else{
        MPI_Request req[2];
        if(params.worm_tail.first >= params.my_start && params.worm_tail.first <= params.my_end){
            first_part = new_coordinates[params.worm_tail.first%params.slices_per_process];
            req[0] = MPI_REQUEST_NULL;
        }
        else
            MPI_Irecv(&first_part[0],params.dimensions,MPI_DOUBLE,params.worm_tail.first/params.slices_per_process,1,local_comm, &req[0]);
        if(params.worm_head.first >= params.my_start && params.worm_head.first <= params.my_end){
            second_part = new_coordinates[params.worm_head.first%params.slices_per_process];
            req[1] = MPI_REQUEST_NULL;
        }
        else
            MPI_Irecv(&second_part[0],params.dimensions,MPI_DOUBLE,params.worm_head.first/params.slices_per_process,2,local_comm, &req[1]);
        MPI_Waitall(2, req, MPI_STATUS_IGNORE);
        std::vector<double> new_action;
        new_action.reserve(params.total_slices);
        new_action.insert(new_action.end(),new_action_p.begin(),new_action_p.end());
        for(int j = 1; j < params.num_workers; ++j){
            MPI_Recv(&new_action_p[0], params.slices_per_process, MPI_DOUBLE,j,0,local_comm,MPI_STATUS_IGNORE);
            new_action.insert(new_action.end(),new_action_p.begin(),new_action_p.end());
        }
        std::vector<double> dist;
        distance(first_part, second_part, dist, params.box_size);
        double r2 = inner_product(dist.begin(),dist.end(), dist.begin(),0.0);
        double rho0 = pow((4* M_PI*params.lambda*disty*params.tau),-params.dimensions/2.)*exp(-1.0/(4.0*params.lambda*disty*params.tau)*r2);
        double na = 0;
        start_slice = (params.worm_tail.first)%params.total_slices;
        for(int j = 0; j < disty; ++j){
            int slice = (start_slice+j)%params.total_slices;
            na += params.tau/2*(new_action[slice]+new_action[(slice+1)%params.total_slices]);
        }
        ac_re = false;
        if(rng.randnormed(1) < (rho0*exp(-na)*exp(params.mu*disty*params.tau))/(params.C*params.Mbar*params.particles*params.total_slices))
            ac_re = true;
    }
}

int Insert::attempt(int &id, Parameters &params, Paths &paths, RNG &rng){
    if(params.worm_on) return 0;
    Moves::attempt(id, params, paths, rng);
    insert_info.push_back(rng.randint(params.total_slices));
    insert_info.push_back(rng.randint(params.Mbar)+1);
    if(params.charged)
        insert_info.push_back(rng.randint(params.charges)*2-1);
    else
        insert_info.push_back(0);
    MPI_Bcast(&insert_info[0], 3, MPI_INT, 0, local_comm);
    new_coordinates.resize(params.slices_per_process);
    new_distances.resize(params.slices_per_process);
    int start_slice = insert_info[0];
    if(start_slice >= params.my_start && start_slice <= params.my_end){
        new_coordinates[start_slice%params.slices_per_process].resize(params.dimensions);
        for(auto &i : new_coordinates[start_slice%params.slices_per_process]){
            if(params.box_size != -1)
                i = rng.randnormed(params.box_size);
            else
                i = rng.randnormed(2)-1;
        }
    }
    for(int i = 1; i <= insert_info[1]; ++i){
        int slice = (start_slice+i)%params.total_slices;
        int slicem1 = (start_slice+i-1)%params.total_slices;
        new_coordinates[slice%params.slices_per_process].resize(params.dimensions);
        if(slice >= params.my_start && slice <= params.my_end){
            if(slicem1 >= params.my_start && slicem1 <= params.my_end)
                new_coordinates[slice%params.slices_per_process] = new_coordinates[slicem1%params.slices_per_process];
            else
                MPI_Recv(&new_coordinates[slice%params.slices_per_process][0], params.dimensions, MPI_DOUBLE, slicem1/params.slices_per_process, i, local_comm, MPI_STATUS_IGNORE);
            double width = sqrt(2*params.lambda*params.tau);
            for(auto &j : new_coordinates[slice%params.slices_per_process])
                j += rng.randgaussian(width);
        }
        else if(slicem1 >= params.my_start && slicem1 <= params.my_end)
                MPI_Send(&new_coordinates[slicem1%params.slices_per_process][0], params.dimensions,MPI_DOUBLE,slice/params.slices_per_process, i, local_comm);
    }
    ++num_attempts;
    check(id, params, paths, rng);
    MPI_Bcast(&ac_re, 1, MPI_INT,0, local_comm);
    if(ac_re){
        ++num_accepts;
        for(int i = 0; i <= insert_info[1]; ++i){
            int slice = (insert_info[0]+i)%params.total_slices;
            if(slice >= params.my_start && slice <= params.my_end){
                put_in_box(new_coordinates[slice%params.slices_per_process], params.box_size);
                paths.worm_advance_tail(params, new_coordinates[slice%params.slices_per_process], new_distances[slice%params.slices_per_process], insert_info[0], insert_info[2]);
            }
            else
                paths.worm_advance_tail(params, std::vector<double>(0), std::vector<std::tuple<int, std::vector<double>, double> >(0), insert_info[0], insert_info[2]);
        }
    }
    insert_info.clear();
    new_coordinates.clear();
    new_distances.clear();
    return 0;
}

void Insert::check(int &id, Parameters &params, Paths &paths, RNG &rng){
    std::vector<double> new_action_p(params.slices_per_process,0);
    std::vector<double> dist(params.dimensions,0);
    std::vector<int> kVec(params.dimensions);
    if(params.potential == 2){
        params.alpha = sqrt(M_PI)*pow(params.particles/params.volume2,1/6.);
        params.coulcut = sqrt(params.p)/params.alpha;
        params.coulcut2 = pow(params.coulcut,2);
        params.kcut = 2.*params.p/params.coulcut;
        params.kcut2 = pow(params.kcut,2);
        params.nmax = floor(params.coulcut/params.box_size);
        params.kmax = ceil(params.kcut/(2.*M_PI/params.box_size));
    }
    int start_slice = insert_info[0];
    for(int i = 0; i <= insert_info[1]; ++i){
        int slice = (start_slice+i)%params.total_slices;
        if(slice >= params.my_start && slice <= params.my_end){
            double new_action_local = 0.0;
            new_distances[slice%params.slices_per_process].reserve(params.particles);
            switch(params.potential){
                case 0:
                    new_action_p[slice%params.slices_per_process] = potential_value_harmonic(inner_product(new_coordinates[slice%params.slices_per_process].begin(), new_coordinates[slice%params.slices_per_process].end(), new_coordinates[slice%params.slices_per_process].begin(), 0.0));
                    break;
                case 1:
#pragma omp parallel for firstprivate(dist) reduction (+:new_action_local)
                    for(int ptcl = 0; ptcl < params.particles; ++ptcl){
                        if(!(ptcl == params.worm_head.second && slice < params.worm_head.first) && !(ptcl == params.worm_tail.second && slice > params.worm_tail.first)){
                            distance(paths.get_coordinate(slice%params.slices_per_process,ptcl),new_coordinates[slice%params.slices_per_process], dist, params.box_size);
                            double r =sqrt(std::inner_product(dist.begin(), dist.end(), dist.begin(), 0.0));
                            new_action_local += potential_value_aziz(r);
#pragma omp critical
                            new_distances[slice%params.slices_per_process].push_back(std::tuple<int,std::vector<double>,double>(paths.get_coordinate_key(slice%params.slices_per_process, ptcl), dist, r));
                        }
                    }
                    new_action_p[slice%params.slices_per_process] += new_action_local;
                    break;
                case 2:
#pragma omp parallel for firstprivate(dist) reduction (+:new_action_local)
                    for(int ptcl = 0; ptcl < params.particles; ++ptcl){
                        if(!(ptcl == params.worm_head.second && slice < params.worm_head.first) && !(ptcl == params.worm_tail.second && slice > params.worm_tail.first)){
                            distance(paths.get_coordinate(slice%params.slices_per_process,ptcl),new_coordinates[slice%params.slices_per_process], dist, params.box_size);
                            double r =sqrt(std::inner_product(dist.begin(), dist.end(), dist.begin(), 0.0));
                            new_action_local += potential_value_coulomb(dist, kVec, paths.charge[ptcl], insert_info[2], params);
#pragma omp critical
                            new_distances[slice%params.slices_per_process].push_back(std::tuple<int,std::vector<double>,double>(paths.get_coordinate_key(slice%params.slices_per_process, ptcl), dist, r));
                        }
                    }
                    dist = std::vector<double>(params.dimensions,0);
                    new_action_p[slice%params.slices_per_process] += new_action_local;
                    new_action_p[slice%params.slices_per_process] += potential_value_coulomb(dist, kVec, insert_info[2], insert_info[2], params);
                    break;
            }
        }
    }
    if(id != 0)
        MPI_Send(&new_action_p[0], params.slices_per_process, MPI_DOUBLE, 0,0,local_comm);
    else{
        std::vector<double> new_action;
        new_action.reserve(params.total_slices);
        new_action.insert(new_action.end(),new_action_p.begin(),new_action_p.end());
        for(int j = 1; j < params.num_workers; ++j){
            MPI_Recv(&new_action_p[0], params.slices_per_process, MPI_DOUBLE,j,0,local_comm,MPI_STATUS_IGNORE);
            new_action.insert(new_action.end(),new_action_p.begin(),new_action_p.end());
        }
        double na = 0;
        for(int j = 0; j < insert_info[1]; ++j){
            int slice = (start_slice+j)%params.total_slices;
            na += params.tau/2*(new_action[slice]+new_action[(slice+1)%params.total_slices]);
        }
        ac_re = false;
        if(rng.randnormed(1) < params.C0*exp(-na)*exp(params.mu*insert_info[1]*params.tau))
            ac_re = true;
    }
}

int Remove::attempt(int &id, Parameters &params, Paths &paths, RNG &rng){
    if(!params.worm_on) return 0;
    if(params.worm_length > params.Mbar + 1) return 1;
    Moves::attempt(id, params, paths, rng);
    check(id, params, paths, rng);
    MPI_Bcast(&ac_re, 1, MPI_INT,0, local_comm);
    ++num_attempts;
    if(ac_re){
        ++num_accepts;
        while(params.worm_on)
            paths.worm_recede_tail(params);
    }
    return 0;
}

void Remove::check(int &id, Parameters &params, Paths &paths, RNG &rng){
    std::vector<double> old_action_p(params.slices_per_process,0);
    std::vector<int> kVec(params.dimensions);
    if(params.potential == 2){
        params.alpha = sqrt(M_PI)*pow(params.particles/params.volume2,1/6.);
        params.coulcut = sqrt(params.p)/params.alpha;
        params.coulcut2 = pow(params.coulcut,2);
        params.kcut = 2.*params.p/params.coulcut;
        params.kcut2 = pow(params.kcut,2);
        params.nmax = floor(params.coulcut/params.box_size);
        params.kmax = ceil(params.kcut/(2.*M_PI/params.box_size));
    }
    int slice = params.worm_head.first;
    int column = params.worm_head.second;
    for(int i = 0; i < params.worm_length; ++i){
        if(slice >= params.my_start && slice <= params.my_end){
            double old_action_local = 0.0;
            switch(params.potential){
                case 0:
                    old_action_p[slice%params.slices_per_process] = potential_value_harmonic(inner_product(paths.get_coordinate(slice%params.slices_per_process,column).begin(),paths.get_coordinate(slice%params.slices_per_process,column).end(),paths.get_coordinate(slice%params.slices_per_process,column).begin(),0.0));
                    break;
                case 1:
#pragma omp parallel for reduction (+:old_action_local)
                    for(int ptcl = 0; ptcl < params.particles; ++ptcl){
                        if(!(paths.broken[ptcl])){
                            if(ptcl != column){
                                old_action_local += potential_value_aziz(paths.get_separation(slice%params.slices_per_process, ptcl, column));
                            }
                        }
                    }
                    old_action_p[slice%params.slices_per_process] += old_action_local;
                    break;
                case 2:
#pragma omp parallel for reduction (+:old_action_local)
                    for(int ptcl = 0; ptcl < params.particles; ++ptcl){
                        if(!(ptcl == params.worm_head.second && slice < params.worm_head.first) && !(ptcl == params.worm_tail.second && slice > params.worm_tail.first)){
                            old_action_local += potential_value_coulomb(paths.get_separation_vector(slice%params.slices_per_process, ptcl, column), kVec, paths.charge[ptcl], paths.charge[column], params);
                        }
                    }
                    old_action_p[slice%params.slices_per_process] += old_action_local;
                    break;
            }
        }
        if(++slice >= params.total_slices){
            slice = 0;
            column = paths.forward_connects[column];
        }
    }
    if(id != 0)
        MPI_Send(&old_action_p[0], params.slices_per_process, MPI_DOUBLE, 0,0,local_comm);
    else{
        ac_re = false;
        std::vector<double> old_action;
        old_action.reserve(params.total_slices);
        old_action.insert(old_action.end(),old_action_p.begin(),old_action_p.end());
        for(int j = 1; j < params.num_workers; ++j){
            MPI_Recv(&old_action_p[0], params.slices_per_process, MPI_DOUBLE,j,0,local_comm,MPI_STATUS_IGNORE);
            old_action.insert(old_action.end(),old_action_p.begin(),old_action_p.end());
        }
        double oa = 0;
        int start_slice = params.worm_head.first;
        if(params.worm_length == 1)
            oa += params.tau*old_action[start_slice];
        else
            for(int j = 0; j < params.worm_length-1; ++j){
                int slice = (start_slice+j)%params.total_slices;
                oa += params.tau/2*(old_action[slice]+old_action[(slice+1)%params.total_slices]);
            }
        if(rng.randnormed(1) < exp(oa)*exp(-params.mu*params.worm_length*params.tau)/params.C0)
            ac_re = true;
    }
}

int Advance_Tail::attempt(int &id, Parameters &params, Paths &paths, RNG &rng){
    if(!params.worm_on) return 0;
    Moves::attempt(id, params, paths, rng);
    M = rng.randint(params.Mbar)+1;
    MPI_Bcast(&M, 1, MPI_INT, 0, local_comm);
    new_coordinates.resize(params.slices_per_process);
    new_distances.resize(params.slices_per_process);
    for(int i = 1; i <= M; ++i){
        int slice = (params.worm_tail.first+i)%params.total_slices;
        int slicem1 = (params.worm_tail.first+i-1)%params.total_slices;
        if(slice >= params.my_start && slice <= params.my_end){
            new_coordinates[slice%params.slices_per_process].resize(params.dimensions);
            if(slicem1 >= params.my_start && slicem1 <= params.my_end){
                if(i == 1)
                    new_coordinates[slicem1%params.slices_per_process] = paths.get_coordinate(slicem1%params.slices_per_process,params.worm_tail.second);
                new_coordinates[slice%params.slices_per_process] = new_coordinates[slicem1%params.slices_per_process];
            }
            else
                MPI_Recv(&new_coordinates[slice%params.slices_per_process][0],params.dimensions,MPI_DOUBLE,slicem1/params.slices_per_process,i, local_comm, MPI_STATUS_IGNORE);
            double width = sqrt(2*params.lambda*params.tau);
            for(auto &j : new_coordinates[slice%params.slices_per_process])
                j += rng.randgaussian(width);
        }
        else if(slicem1 >= params.my_start && slicem1 <= params.my_end){
            if(i == 1)
                new_coordinates[slicem1%params.slices_per_process] = paths.get_coordinate(slicem1%params.slices_per_process,params.worm_tail.second);
            MPI_Send(&new_coordinates[slicem1%params.slices_per_process][0],params.dimensions,MPI_DOUBLE,slice/params.slices_per_process,i, local_comm);
        }
    }
    check(id, params, paths, rng);
    MPI_Bcast(&ac_re, 1, MPI_INT,0, local_comm);
    ++num_attempts;
    if(ac_re){
        ++num_accepts;
        int start_slice = params.worm_tail.first;
        for(int i = 1; i <= M; ++i){
            int slice = (start_slice+i)%params.total_slices;
            if(slice >= params.my_start && slice <= params.my_end){
                put_in_box(new_coordinates[slice%params.slices_per_process], params.box_size);
                paths.worm_advance_tail(params, new_coordinates[slice%params.slices_per_process],new_distances[slice%params.slices_per_process]);
            }
            else
                paths.worm_advance_tail(params, std::vector<double>(0), std::vector<std::tuple<int, std::vector<double>, double> >(0));
        }
    }
    new_coordinates.clear();
    new_distances.clear();
    return 0;
}

void Advance_Tail::check(int &id, Parameters &params, Paths &paths, RNG &rng){
    std::vector<double> new_action_p(params.slices_per_process,0);
    std::vector<double> dist(params.dimensions,0);
    std::vector<int> kVec(params.dimensions);
    if(params.potential == 2){
        params.alpha = sqrt(M_PI)*pow(params.particles/params.volume2,1/6.);
        params.coulcut = sqrt(params.p)/params.alpha;
        params.coulcut2 = pow(params.coulcut,2);
        params.kcut = 2.*params.p/params.coulcut;
        params.kcut2 = pow(params.kcut,2);
        params.nmax = floor(params.coulcut/params.box_size);
        params.kmax = ceil(params.kcut/(2.*M_PI/params.box_size));
    }
    int start_slice = params.worm_tail.first;
    for(int i = 0; i <= M; ++i){
        int slice = (start_slice+i)%params.total_slices;
        if(slice >= params.my_start && slice <= params.my_end){
            double new_action_local = 0;
            new_distances[slice%params.slices_per_process].reserve(params.particles);
            switch(params.potential){
                case 0:
                    new_action_p[slice%params.slices_per_process] = potential_value_harmonic(inner_product(new_coordinates[slice%params.slices_per_process].begin(), new_coordinates[slice%params.slices_per_process].end(), new_coordinates[slice%params.slices_per_process].begin(), 0.0));
                    break;
                case 1:
#pragma omp parallel for firstprivate(dist) reduction (+:new_action_local)
                    for(int ptcl = 0; ptcl < params.particles; ++ptcl){
                        if(!(ptcl == params.worm_head.second && slice < params.worm_head.first) && !(ptcl == params.worm_tail.second && slice >= params.worm_tail.first)){
                            distance(paths.get_coordinate(slice%params.slices_per_process,ptcl), new_coordinates[slice%params.slices_per_process], dist, params.box_size);
                            double r = sqrt(std::inner_product(dist.begin(), dist.end(), dist.begin(), 0.0));
                            new_action_local  += potential_value_aziz(r);
#pragma omp critical
                            new_distances[slice%params.slices_per_process].push_back(std::tuple<int,std::vector<double>,double>(paths.get_coordinate_key(slice%params.slices_per_process, ptcl), dist, r));
                        }
                    }
                    new_action_p[slice%params.slices_per_process] += new_action_local;
                    break;
                case 2:
#pragma omp parallel for firstprivate(dist) reduction (+:new_action_local)
                    for(int ptcl = 0; ptcl < params.particles; ++ptcl){
                        if(!(ptcl == params.worm_head.second && slice < params.worm_head.first) && !(ptcl == params.worm_tail.second && slice >= params.worm_tail.first)){
                            distance(paths.get_coordinate(slice%params.slices_per_process,ptcl), new_coordinates[slice%params.slices_per_process], dist, params.box_size);
                            double r = sqrt(std::inner_product(dist.begin(), dist.end(), dist.begin(), 0.0));
                            new_action_local += potential_value_coulomb(dist, kVec, paths.charge[ptcl], paths.charge[params.worm_tail.second], params);
#pragma omp critical
                            new_distances[slice%params.slices_per_process].push_back(std::tuple<int,std::vector<double>,double>(paths.get_coordinate_key(slice%params.slices_per_process, ptcl), dist, r));
                        }
                    }
                    new_action_p[slice%params.slices_per_process] += new_action_local;
                    dist = std::vector<double>(params.dimensions,0);
                    new_action_p[slice%params.slices_per_process] += potential_value_coulomb(dist, kVec, paths.charge[params.worm_tail.second], paths.charge[params.worm_tail.second], params);
                    break;
            }
        }
    }
    if(id != 0)
        MPI_Send(&new_action_p[0], params.slices_per_process, MPI_DOUBLE, 0,0,local_comm);
    else{
        ac_re = false;
        std::vector<double> new_action;
        new_action.reserve(params.total_slices);
        new_action.insert(new_action.end(),new_action_p.begin(),new_action_p.end());
        for(int j = 1; j < params.num_workers; ++j){
            MPI_Recv(&new_action_p[0], params.slices_per_process, MPI_DOUBLE,j,0,local_comm,MPI_STATUS_IGNORE);
            new_action.insert(new_action.end(),new_action_p.begin(),new_action_p.end());
        }
        double na = 0;
        for(int j = 0; j < M; ++j){
            int slice = (start_slice+j)%params.total_slices;
            na += params.tau/2*(new_action[slice]+new_action[(slice+1)%params.total_slices]);
        }
        if(rng.randnormed(1) < exp(-na)*exp(params.mu*M*params.tau))
            ac_re = true;
    }
}

int Advance_Head::attempt(int &id, Parameters &params, Paths &paths, RNG &rng){
    if(!params.worm_on) return 0;
    Moves::attempt(id, params, paths, rng);
    M = rng.randint(params.Mbar)+1;
    MPI_Bcast(&M, 1, MPI_INT, 0, local_comm);
    new_coordinates.resize(params.slices_per_process);
    new_distances.resize(params.slices_per_process);
    for(int i = 1; i <= M; ++i){
        int slice = (params.worm_head.first-i+params.total_slices)%params.total_slices;
        int slicep1 = (params.worm_head.first-i+1+params.total_slices)%params.total_slices;
        if(slice >= params.my_start && slice <= params.my_end){
            new_coordinates[slice%params.slices_per_process].resize(params.dimensions);
            if(slicep1 >= params.my_start && slicep1 <= params.my_end){
                if(i == 1)
                    new_coordinates[slicep1%params.slices_per_process] = paths.get_coordinate(slicep1%params.slices_per_process,params.worm_head.second);
                new_coordinates[slice%params.slices_per_process] = new_coordinates[slicep1%params.slices_per_process];
            }
            else
                MPI_Recv(&new_coordinates[slice%params.slices_per_process][0],params.dimensions,MPI_DOUBLE,slicep1/params.slices_per_process,i, local_comm, MPI_STATUS_IGNORE);
            double width = sqrt(2*params.lambda*params.tau);
            for(auto &j : new_coordinates[slice%params.slices_per_process])
                j += rng.randgaussian(width);
        }
        else if(slicep1 >= params.my_start && slicep1 <= params.my_end){
            if(i == 1)
                new_coordinates[slicep1%params.slices_per_process] = paths.get_coordinate(slicep1%params.slices_per_process,params.worm_head.second);
            MPI_Send(&new_coordinates[slicep1%params.slices_per_process][0],params.dimensions,MPI_DOUBLE,slice/params.slices_per_process,i, local_comm);
        }
    }
    check(id, params, paths, rng);
    MPI_Bcast(&ac_re, 1, MPI_INT,0, local_comm);
    ++num_attempts;
    if(ac_re){
        ++num_accepts;
        int start_slice = params.worm_head.first;
        for(int i = 1; i <= M; ++i){
            int slice = (start_slice-i+params.total_slices)%params.total_slices;
            if(slice >= params.my_start && slice <= params.my_end){
                put_in_box(new_coordinates[slice%params.slices_per_process], params.box_size);
                paths.worm_advance_head(params, new_coordinates[slice%params.slices_per_process],new_distances[slice%params.slices_per_process]);
            }
            else
                paths.worm_advance_head(params, std::vector<double>(0),std::vector<std::tuple<int, std::vector<double>, double> >(0));
        }
    }
    new_coordinates.clear();
    new_distances.clear();
    return 0;
}

void Advance_Head::check(int &id, Parameters &params, Paths &paths, RNG &rng){
    int start_slice = params.worm_head.first;
    std::vector<double> new_action_p(params.slices_per_process,0);
    std::vector<double> dist(params.dimensions,0);
    std::vector<int> kVec(params.dimensions);
    if(params.potential == 2){
        params.alpha = sqrt(M_PI)*pow(params.particles/params.volume2,1/6.);
        params.coulcut = sqrt(params.p)/params.alpha;
        params.coulcut2 = pow(params.coulcut,2);
        params.kcut = 2.*params.p/params.coulcut;
        params.kcut2 = pow(params.kcut,2);
        params.nmax = floor(params.coulcut/params.box_size);
        params.kmax = ceil(params.kcut/(2.*M_PI/params.box_size));
    }
    for(int i = 0; i <= M; ++i){
        int slice = (start_slice-i+params.total_slices)%params.total_slices;
        if(slice >= params.my_start && slice <= params.my_end){
            double new_action_local = 0;
            new_distances[slice%params.slices_per_process].reserve(params.particles);
            switch(params.potential){
                case 0:
                    new_action_p[slice%params.slices_per_process] = potential_value_harmonic(inner_product(new_coordinates[slice%params.slices_per_process].begin(), new_coordinates[slice%params.slices_per_process].end(), new_coordinates[slice%params.slices_per_process].begin(), 0.0));
                    break;
                case 1:
#pragma omp parallel for firstprivate(dist) reduction (+:new_action_local)
                    for(int ptcl = 0; ptcl < params.particles; ++ptcl){
                        if(!(ptcl == params.worm_head.second && slice <= params.worm_head.first) && !(ptcl == params.worm_tail.second && slice > params.worm_tail.first)){
                            distance(paths.get_coordinate(slice%params.slices_per_process,ptcl), new_coordinates[slice%params.slices_per_process], dist, params.box_size);
                            double r = sqrt(std::inner_product(dist.begin(), dist.end(), dist.begin(), 0.0));
                            new_action_local += potential_value_aziz(r);
#pragma omp critical
                            new_distances[slice%params.slices_per_process].push_back(std::tuple<int,std::vector<double>,double>(paths.get_coordinate_key(slice%params.slices_per_process, ptcl), dist, r));
                        }
                    }
                    new_action_p[slice%params.slices_per_process] += new_action_local;
                    break;
                case 2:
#pragma omp parallel for firstprivate(dist) reduction (+:new_action_local)
                    for(int ptcl = 0; ptcl < params.particles; ++ptcl){
                        if(!(ptcl == params.worm_head.second && slice <= params.worm_head.first) && !(ptcl == params.worm_tail.second && slice > params.worm_tail.first)){
                            distance(paths.get_coordinate(slice%params.slices_per_process,ptcl), new_coordinates[slice%params.slices_per_process], dist, params.box_size);
                            double r = sqrt(std::inner_product(dist.begin(), dist.end(), dist.begin(), 0.0));
                            new_action_local += potential_value_coulomb(dist, kVec, paths.charge[ptcl], paths.charge[params.worm_head.second], params);
#pragma omp critical
                            new_distances[slice%params.slices_per_process].push_back(std::tuple<int,std::vector<double>,double>(paths.get_coordinate_key(slice%params.slices_per_process, ptcl), dist, r));
                        }
                    }
                    new_action_p[slice%params.slices_per_process] += new_action_local;
                    dist = std::vector<double>(params.dimensions,0);
                    new_action_p[slice%params.slices_per_process] += potential_value_coulomb(dist, kVec, paths.charge[params.worm_head.second], paths.charge[params.worm_head.second], params);
                    break;
            }
        }
    }
    if(id != 0)
        MPI_Send(&new_action_p[0], params.slices_per_process, MPI_DOUBLE, 0,0,local_comm);
    else{
        ac_re = false;
        std::vector<double> new_action;
        new_action.reserve(params.total_slices);
        new_action.insert(new_action.end(),new_action_p.begin(),new_action_p.end());
        for(int j = 1; j < params.num_workers; ++j){
            MPI_Recv(&new_action_p[0], params.slices_per_process, MPI_DOUBLE,j,0,local_comm,MPI_STATUS_IGNORE);
            new_action.insert(new_action.end(),new_action_p.begin(),new_action_p.end());
        }
        double na = 0;
        for(int j = 0; j < M; ++j){
            int slice = (start_slice-j+params.total_slices)%params.total_slices;
            na += params.tau/2*(new_action[slice]+new_action[(slice-1+params.total_slices)%params.total_slices]);
        }
        if(rng.randnormed(1) < exp(-na)*exp(params.mu*M*params.tau))
            ac_re = true;
    }
}

int Recede_Head::attempt(int &id, Parameters &params, Paths &paths, RNG &rng){
    if(!params.worm_on) return 0;
    Moves::attempt(id, params, paths, rng);
    M = rng.randint(params.Mbar)+1;
    if(M >= params.worm_length)
        M = params.worm_length;
    MPI_Bcast(&M, 1, MPI_INT, 0, local_comm);
    check(id, params, paths, rng);
    MPI_Bcast(&ac_re, 1, MPI_INT,0, local_comm);
    ++num_attempts;
    if(ac_re){
        ++num_accepts;
        for(int i = 0; i < M; ++i)
            paths.worm_recede_head(params);
    }
    return 0;
}

void Recede_Head::check(int &id, Parameters &params, Paths &paths, RNG &rng){
    std::vector<double> old_action_p(params.slices_per_process,0);
    std::vector<int> kVec(params.dimensions);
    if(params.potential == 2){
        params.alpha = sqrt(M_PI)*pow(params.particles/params.volume2,1/6.);
        params.coulcut = sqrt(params.p)/params.alpha;
        params.coulcut2 = pow(params.coulcut,2);
        params.kcut = 2.*params.p/params.coulcut;
        params.kcut2 = pow(params.kcut,2);
        params.nmax = floor(params.coulcut/params.box_size);
        params.kmax = ceil(params.kcut/(2.*M_PI/params.box_size));
    }
    int start_slice = params.worm_head.first;
    int col = params.worm_head.second;
    for(int i = 0; i < M; ++i){
        int slice = (start_slice+i)%params.total_slices;
        if(start_slice+i == params.total_slices)
            col = paths.forward_connects[col];
        if(slice >= params.my_start && slice <= params.my_end){
            double old_action_local = 0;
            switch(params.potential){
                case 0:
                    old_action_p[slice%params.slices_per_process] = potential_value_harmonic(inner_product(paths.get_coordinate(slice%params.slices_per_process,col).begin(),paths.get_coordinate(slice%params.slices_per_process,col).end(),paths.get_coordinate(slice%params.slices_per_process,col).begin(),0.0));
                    break;
                case 1:
#pragma omp parallel for reduction (+:old_action_local)
                    for(int ptcl = 0; ptcl < params.particles; ++ptcl){
                        if(!(ptcl == params.worm_head.second && slice < params.worm_head.first) && !(ptcl == params.worm_tail.second && slice > params.worm_tail.first)){
                            if(ptcl != col){
                                old_action_local += potential_value_aziz(paths.get_separation(slice%params.slices_per_process, ptcl, col));
                            }
                        }
                    }
                    old_action_p[slice%params.slices_per_process] += old_action_local;
                    break;
                case 2:
#pragma omp parallel for reduction (+:old_action_local)
                    for(int ptcl = 0; ptcl < params.particles; ++ptcl){
                        if(!(ptcl == params.worm_head.second && slice < params.worm_head.first) && !(ptcl == params.worm_tail.second && slice > params.worm_tail.first)){
                            old_action_local += potential_value_coulomb(paths.get_separation_vector(slice%params.slices_per_process, ptcl, col), kVec, paths.charge[ptcl], paths.charge[col], params);
                        }
                    }
                    old_action_p[slice%params.slices_per_process] += old_action_local;
                    break;
            }
        }
    }
    if(id != 0)
        MPI_Send(&old_action_p[0], params.slices_per_process, MPI_DOUBLE, 0,0,local_comm);
    else{
        ac_re = false;
        std::vector<double> old_action;
        old_action.reserve(params.total_slices);
        for(int j = 1; j < params.num_workers; ++j){
            MPI_Recv(&old_action_p[0], params.slices_per_process, MPI_DOUBLE,j,0,local_comm,MPI_STATUS_IGNORE);
            old_action.insert(old_action.end(),old_action_p.begin(),old_action_p.end());
        }
        double oa = 0;
        for(int j = 0; j < M-1; ++j){
            int slice = (start_slice+j)%params.total_slices;
            oa += params.tau/2*(old_action[slice]+old_action[(slice+1)%params.total_slices]);
        }
        double rn = rng.randnormed(1);
        if(rn < exp(oa)*exp(-params.mu*M*params.tau))
            ac_re = true;
    }
}

int Recede_Tail::attempt(int &id, Parameters &params, Paths &paths, RNG &rng){
    if(!params.worm_on) return 0;
    Moves::attempt(id, params, paths, rng);
    M = rng.randint(params.Mbar)+1;
    if(M > params.worm_length) M = params.worm_length;
    MPI_Bcast(&M, 1, MPI_INT, 0, local_comm);
    check(id, params, paths, rng);
    MPI_Bcast(&ac_re, 1, MPI_INT,0, local_comm);
    ++num_attempts;
    if(ac_re){
        ++num_accepts;
        for(int i = 0; i < M; ++i)
            paths.worm_recede_tail(params);
    }
    return 0;
}

void Recede_Tail::check(int &id, Parameters &params, Paths &paths, RNG &rng){
    std::vector<double> old_action_p(params.slices_per_process,0);
    std::vector<int> kVec(params.dimensions);
    if(params.potential == 2){
        params.alpha = sqrt(M_PI)*pow(params.particles/params.volume2,1/6.);
        params.coulcut = sqrt(params.p)/params.alpha;
        params.coulcut2 = pow(params.coulcut,2);
        params.kcut = 2.*params.p/params.coulcut;
        params.kcut2 = pow(params.kcut,2);
        params.nmax = floor(params.coulcut/params.box_size);
        params.kmax = ceil(params.kcut/(2.*M_PI/params.box_size));
    }
    int start_slice = params.worm_tail.first;
    int col = params.worm_tail.second;
    for(int i = 0; i < M; i++){
        int slice = (start_slice-i+params.total_slices)%params.total_slices;
        if(start_slice - i == -1)
            col = paths.backward_connects[col];
        if(slice >= params.my_start && slice <= params.my_end){
            double old_action_local = 0;
            switch(params.potential){
                case 0:
                    old_action_p[slice%params.slices_per_process] = potential_value_harmonic(inner_product(paths.get_coordinate(slice%params.slices_per_process,col).begin(),paths.get_coordinate(slice%params.slices_per_process,col).end(),paths.get_coordinate(slice%params.slices_per_process,col).begin(),0.0));
                    break;
                case 1:
#pragma omp parallel for reduction (+:old_action_local)
                    for(int ptcl = 0; ptcl < params.particles; ++ptcl){
                        if(!(ptcl == params.worm_head.second && slice < params.worm_head.first) && !(ptcl == params.worm_tail.second && slice > params.worm_tail.first)){
                            if(ptcl != col){
                                old_action_local += potential_value_aziz(paths.get_separation(slice%params.slices_per_process, ptcl, col));
                            }
                        }
                    }
                    old_action_p[slice%params.slices_per_process] += old_action_local;
                    break;
                case 2:
#pragma omp parallel for reduction (+:old_action_local)
                    for(int ptcl = 0; ptcl < params.particles; ++ptcl){
                        if(!(ptcl == params.worm_head.second && slice < params.worm_head.first) && !(ptcl == params.worm_tail.second && slice > params.worm_tail.first)){
                            old_action_local += potential_value_coulomb(paths.get_separation_vector(slice%params.slices_per_process, ptcl, col), kVec, paths.charge[ptcl], paths.charge[col], params);
                        }
                    }
                    old_action_p[slice%params.slices_per_process] += old_action_local;
                    break;
            }
        }
    }
    if(id != 0)
        MPI_Send(&old_action_p[0], params.slices_per_process, MPI_DOUBLE, 0,0,local_comm);
    else{
        ac_re = false;
        std::vector<double> old_action;
        old_action.reserve(params.total_slices);
        old_action.insert(old_action.end(),old_action_p.begin(),old_action_p.end());
        for(int j = 1; j < params.num_workers; ++j){
            MPI_Recv(&old_action_p[0], params.slices_per_process, MPI_DOUBLE,j,0,local_comm,MPI_STATUS_IGNORE);
            old_action.insert(old_action.end(),old_action_p.begin(),old_action_p.end());
        }
        double oa = 0;
        for(int j = 0; j < M-1; ++j){
            int slice = (start_slice-j+params.total_slices)%params.total_slices;
            oa += params.tau/2*(old_action[slice]+old_action[(slice-1+params.total_slices)%params.total_slices]);
        }
        if(rng.randnormed(1) < exp(oa)*exp(-params.mu*M*params.tau))
            ac_re = true;
    }
}

int Swap_Tail::attempt(int &id, Parameters &params, Paths &paths, RNG &rng){
    if(!params.worm_on) return 0;
    Moves::attempt(id, params, paths, rng);
    int start_slice = params.worm_tail.first;
    int swap_slice = (start_slice + params.Mbar) % params.total_slices;
    keep_going = true;
    sigma_I = 0;
    sigma_Z = 0;
    std::vector<double> start(params.dimensions,0);
    std::vector<double> end(params.dimensions,0);
    if(start_slice >= params.my_start && start_slice <= params.my_end){
        int col = params.worm_tail.second;
        int sslice = start_slice%params.slices_per_process;
        start = paths.get_coordinate(sslice, col);
        if(swap_slice >= params.my_start && swap_slice <= params.my_end){
            int eslice = swap_slice%params.slices_per_process;
            std::vector<int> possible_swaps = paths.get_nearest_neighbors(eslice, start);
            if(!possible_swaps.empty()){
                for(int poss = possible_swaps.size() - 1; poss >= 0; --poss){
                    if(paths.charge[possible_swaps[poss]] != paths.charge[params.worm_tail.second])
                        possible_swaps.erase(possible_swaps.begin() + poss);
                }
            }
            if(possible_swaps.size() == 0)
                keep_going = false;
            else{
                std::vector<double> rho0s(possible_swaps.size(),0);
                std::vector<double> dist(params.dimensions,0);
                for(auto &i : possible_swaps){
                    distance(start, paths.get_coordinate(eslice, i), dist, params.box_size);
                    double r2 = inner_product(dist.begin(), dist.end(), dist.begin(),0.0);
                    rho0s[&i - &possible_swaps[0]] = pow((4* M_PI*params.lambda*params.Mbar*params.tau),-params.dimensions/2.)*exp(-1.0/(4.0*params.lambda*params.Mbar*params.tau)*r2);
                }
                for (auto &r : rho0s)
                    sigma_I += r;
                double total = 0;
                double rn = rng.randnormed(sigma_I);
                for(auto &r: rho0s){
                    total += r;
                    if(rn <= total){
                        choice = possible_swaps[&r - &rho0s[0]];
                        break;
                    }
                }
                if(paths.broken[choice]){
                    int back_part = choice;
                    if(start_slice + params.Mbar >= params.total_slices)
                        back_part = paths.backward_connects[choice];
                    if(back_part == -1 || (back_part == params.worm_head.second && start_slice < params.worm_head.first))
                        keep_going = false;
                }
                if(keep_going){
                    int choice_back = choice;
                    if(start_slice + params.Mbar >= params.total_slices)
                        choice_back = paths.backward_connects[choice];
                    std::vector<double> cb_loc(params.dimensions,0);
                    end = paths.get_coordinate(eslice, choice);
                    cb_loc = paths.get_coordinate(sslice,choice_back);
                    if((keep_going = paths.are_neighbors(cb_loc, end))){
                        std::vector<int> possible_swaps_2 = paths.get_nearest_neighbors(eslice, cb_loc);
                        if(!possible_swaps_2.empty()){
                            for(int poss = possible_swaps_2.size() - 1; poss >= 0; --poss){
                                if(paths.charge[possible_swaps_2[poss]] != paths.charge[params.worm_tail.second])
                                    possible_swaps_2.erase(possible_swaps_2.begin() + poss);
                            }
                        }
                        std::vector<double> rho0s_2(possible_swaps_2.size(),0);
                        for(auto &i : possible_swaps_2){
                            distance(cb_loc, paths.get_coordinate(eslice, i), dist, params.box_size);
                            double r2 = inner_product(dist.begin(),dist.end(), dist.begin(),0.0);
                            rho0s_2[&i - &possible_swaps_2[0]] = pow((4* M_PI*params.lambda*params.Mbar*params.tau),-params.dimensions/2.)*exp(-1.0/(4.0*params.lambda*params.Mbar*params.tau)*r2);
                        }
                        for (auto &r : rho0s_2)
                            sigma_Z += r;
                    }
                }
            }
        }
        else{
            MPI_Send(&start[0], params.dimensions, MPI_DOUBLE,swap_slice/params.slices_per_process, 0, local_comm);
            MPI_Recv(&keep_going, 1, MPI_INT, swap_slice/params.slices_per_process, 1, local_comm, MPI_STATUS_IGNORE);
            if(keep_going){
                MPI_Recv(&keep_going, 1, MPI_INT, swap_slice/params.slices_per_process, 2, local_comm, MPI_STATUS_IGNORE);
                if(keep_going){
                    MPI_Recv(&choice, 1, MPI_INT, swap_slice/params.slices_per_process, 3, local_comm, MPI_STATUS_IGNORE);
                    int choice_back = choice;
                    if(start_slice + params.Mbar >= params.total_slices){
                        choice_back = paths.backward_connects[choice];
                    }
                    std::vector<double> cb_loc(params.dimensions,0);
                    cb_loc = paths.get_coordinate(sslice,choice_back);
                    MPI_Send(&cb_loc[0], params.dimensions, MPI_DOUBLE,swap_slice/params.slices_per_process, 4, local_comm);
                }
            }
        }
    }
    else if(swap_slice >= params.my_start && swap_slice <= params.my_end){
        std::vector<double> tail(params.dimensions,0);
        MPI_Recv(&tail[0], params.dimensions, MPI_DOUBLE, start_slice/params.slices_per_process, 0, local_comm, MPI_STATUS_IGNORE);
        int slice = swap_slice%params.slices_per_process;
        std::vector<int> possible_swaps = paths.get_nearest_neighbors(slice, tail);
        if(!possible_swaps.empty()){
            for(int poss = possible_swaps.size() - 1; poss >= 0; --poss){
                if(paths.charge[possible_swaps[poss]] != paths.charge[params.worm_tail.second])
                    possible_swaps.erase(possible_swaps.begin() + poss);
            }
        }
        if(possible_swaps.size() == 0)
            keep_going = false;
        MPI_Send(&keep_going, 1, MPI_INT, start_slice/params.slices_per_process, 1, local_comm);
        if(keep_going){
            std::vector<double> rho0s(possible_swaps.size(),0);
            std::vector<double> dist(params.dimensions,0);
            for(auto &i : possible_swaps){
                distance(tail, paths.get_coordinate(slice, i), dist, params.box_size);
                double r2 = inner_product(dist.begin(),dist.end(), dist.begin(),0.0);
                rho0s[&i - &possible_swaps[0]] = pow((4* M_PI*params.lambda*params.Mbar*params.tau),-params.dimensions/2.)*exp(-1.0/(4.0*params.lambda*params.Mbar*params.tau)*r2);
            }
            for (auto &r : rho0s)
                sigma_I += r;
            double total = 0;
            double rn = rng.randnormed(sigma_I);
            for(auto &r : rho0s){
                total += r;
                if(rn <= total){
                    choice = possible_swaps[&r - &rho0s[0]];
                    break;
                }
            }
            if(paths.broken[choice]){
                int back_part = choice;
                if(start_slice + params.Mbar >= params.total_slices)
                    back_part = paths.backward_connects[choice];
                if(back_part == -1 || (back_part == params.worm_head.second && start_slice < params.worm_head.first))
                    keep_going = false;
            }
            MPI_Send(&keep_going, 1, MPI_INT, start_slice/params.slices_per_process, 2, local_comm);
            if(keep_going){
                MPI_Send(&choice, 1, MPI_INT, start_slice/params.slices_per_process, 3, local_comm);
                end = paths.get_coordinate(slice, choice);
                std::vector<double> cb(params.dimensions);
                MPI_Recv(&cb[0], params.dimensions, MPI_DOUBLE, start_slice/params.slices_per_process, 4, local_comm, MPI_STATUS_IGNORE);
                if((keep_going = paths.are_neighbors(cb, end))){
                    std::vector<int> possible_swaps_2 = paths.get_nearest_neighbors(slice, cb);
                    if(!possible_swaps_2.empty()){
                        for(int poss = possible_swaps_2.size() - 1; poss >= 0; --poss){
                            if(paths.charge[possible_swaps_2[poss]] != paths.charge[params.worm_tail.second])
                                possible_swaps_2.erase(possible_swaps_2.begin() + poss);
                        }
                    }
                    std::vector<double> rho0s_2(possible_swaps_2.size(),0);
                    for(auto &i : possible_swaps_2){
                        distance(end, paths.get_coordinate(slice, i), dist, params.box_size);
                        double r2 = inner_product(dist.begin(),dist.end(), dist.begin(),0.0);
                        rho0s_2[&i - &possible_swaps_2[0]] = pow((4* M_PI*params.lambda*params.Mbar*params.tau),-params.dimensions/2.)*exp(-1.0/(4.0*params.lambda*params.Mbar*params.tau)*r2);
                    }
                    for (auto &r : rho0s_2)
                        sigma_Z += r;
                }
            }
        }
    }
    MPI_Bcast(&keep_going, 1, MPI_INT, swap_slice/params.slices_per_process, local_comm);
    if(!keep_going) return 1;
    MPI_Bcast(&choice, 1, MPI_INT, swap_slice/params.slices_per_process, local_comm);
    MPI_Bcast(&sigma_I, 1, MPI_DOUBLE, swap_slice/params.slices_per_process, local_comm);
    MPI_Bcast(&sigma_Z, 1, MPI_DOUBLE, swap_slice/params.slices_per_process, local_comm);
    MPI_Bcast(&end[0], params.dimensions, MPI_DOUBLE, swap_slice/params.slices_per_process, local_comm);
    new_coordinates.resize(params.slices_per_process);
    new_distances.reserve(params.Mbar*params.particles);
    for(int i = 1; i < params.Mbar; ++i){
        int slice = (start_slice+i)%params.total_slices;
        int slicem1 = (start_slice+i-1+params.total_slices)%params.total_slices;
        if(swap_slice >= params.my_start && swap_slice <= params.my_end){
            new_coordinates[swap_slice%params.slices_per_process] = end;
        }
        if(slice >= params.my_start && slice <= params.my_end){
            if(slicem1 >= params.my_start && slicem1 <= params.my_end){
                if(i == 1)
                    new_coordinates[slicem1%params.slices_per_process] = start;
                start = new_coordinates[slicem1%params.slices_per_process];
            }
            else
                MPI_Recv(&start[0], params.dimensions, MPI_DOUBLE, slicem1/params.slices_per_process, i, local_comm, MPI_STATUS_IGNORE);
            average_loc_weighted(start, end, new_coordinates[slice%params.slices_per_process], params.box_size, 1, params.Mbar-i);
            double tau_w = params.tau*params.Mbar-i;
            double width = sqrt(2.*params.lambda/(1./params.tau+1./tau_w));
            for(auto &j : new_coordinates[slice%params.slices_per_process])
                j += rng.randgaussian(width);
            put_in_box(new_coordinates[slice%params.slices_per_process], params.box_size);
        }
        else{
            if(slicem1 >= params.my_start && slicem1 <= params.my_end){
                if(i == 1)
                    new_coordinates[slicem1%params.slices_per_process] = start;
                MPI_Send(&new_coordinates[slicem1%params.slices_per_process][0], params.dimensions,MPI_DOUBLE,slice/params.slices_per_process, i, local_comm);
            }
        }
    }
    check(id, params, paths, rng);
    MPI_Bcast(&ac_re, 1, MPI_INT,0, local_comm);
    ++num_attempts;
    if(ac_re){
        ++num_accepts;
        paths.swap_into_tail(params, choice, params.Mbar, new_coordinates);
        paths.update_separations(new_distances);
    }
    new_coordinates.clear();
    new_distances.clear();
    return 0;
}

void Swap_Tail::check(int &id, Parameters &params, Paths &paths, RNG &rng){
    std::vector<double> old_action_p(params.slices_per_process,0);
    std::vector<double> new_action_p(params.slices_per_process,0);
    std::vector<int> kVec(params.dimensions);
    std::vector<double> dist(params.dimensions);
    if(params.potential == 2){
        params.alpha = sqrt(M_PI)*pow(params.particles/params.volume2,1/6.);
        params.coulcut = sqrt(params.p)/params.alpha;
        params.coulcut2 = pow(params.coulcut,2);
        params.kcut = 2.*params.p/params.coulcut;
        params.kcut2 = pow(params.kcut,2);
        params.nmax = floor(params.coulcut/params.box_size);
        params.kmax = ceil(params.kcut/(2.*M_PI/params.box_size));
    }
    int start_slice = params.worm_tail.first;
    for(int i = 0; i <= params.Mbar; ++i){
        int slice = (start_slice+i)%params.total_slices;
        int ptcl1 = choice;
        if(slice >= params.my_start && slice <= params.my_end){
            if(start_slice+params.Mbar >= params.total_slices && start_slice+i < params.total_slices)
                ptcl1 = paths.backward_connects[ptcl1];
            double old_action_local = 0;
            double new_action_local = 0;
            switch(params.potential){
                case 0:
                {
                    new_action_p[slice%params.slices_per_process] = potential_value_harmonic(inner_product(new_coordinates[slice%params.slices_per_process].begin(), new_coordinates[slice%params.slices_per_process].end(), new_coordinates[slice%params.slices_per_process].begin(),0.0));
                    old_action_p[slice%params.slices_per_process] = potential_value_harmonic(inner_product(paths.get_coordinate(slice%params.slices_per_process, ptcl1).begin(), paths.get_coordinate(slice%params.slices_per_process, ptcl1).end(), paths.get_coordinate(slice%params.slices_per_process, ptcl1).begin(), 0.0));
                    break;
                }
                case 1:
#pragma omp parallel for firstprivate(dist) reduction(+:old_action_local,new_action_local)
                    for(int ptcl = 0; ptcl < params.particles; ++ptcl){
                        if(!(ptcl == params.worm_head.second && slice < params.worm_head.first) && !(ptcl == params.worm_tail.second && slice > params.worm_tail.first)){
                            if(i != 0 && ptcl != ptcl1){
                                old_action_local += potential_value_aziz(paths.get_separation(slice%params.slices_per_process, ptcl, ptcl1));
                                distance(paths.get_coordinate(slice%params.slices_per_process,ptcl), new_coordinates[slice%params.slices_per_process], dist, params.box_size);
                                double r = sqrt(std::inner_product(dist.begin(), dist.end(), dist.begin(), 0.0));
                                new_action_local += potential_value_aziz(r);
#pragma omp critical
                                new_distances.push_back(std::tuple<std::pair<int,int>,std::vector<double>,double>(std::pair<int,int>(paths.get_coordinate_key(slice%params.slices_per_process, ptcl),paths.get_coordinate_key(slice%params.slices_per_process, ptcl1)),  dist, r));
                            }
                            else if(i == 0){
                                if(ptcl != ptcl1)
                                    old_action_local += potential_value_aziz(paths.get_separation(slice%params.slices_per_process, ptcl, ptcl1));
                                distance(paths.get_coordinate(slice%params.slices_per_process,ptcl), new_coordinates[slice%params.slices_per_process], dist, params.box_size);
                                double r = sqrt(std::inner_product(dist.begin(), dist.end(), dist.begin(), 0.0));
                                if(ptcl != params.worm_tail.second)
                                    new_action_local += potential_value_aziz(r);
                            }
                        }
                    }
                    old_action_p[slice%params.slices_per_process] += old_action_local;
                    new_action_p[slice%params.slices_per_process] += new_action_local;
                    break;
                case 2:
#pragma omp parallel for firstprivate(dist) reduction(+:old_action_local,new_action_local)
                    for(int ptcl = 0; ptcl < params.particles; ++ptcl){
                        if(!(ptcl == params.worm_head.second && slice < params.worm_head.first) && !(ptcl == params.worm_tail.second && slice > params.worm_tail.first)){
                            old_action_local += potential_value_coulomb(paths.get_separation_vector(slice%params.slices_per_process, ptcl, ptcl1), kVec, paths.charge[ptcl], paths.charge[ptcl1], params);
                            distance(paths.get_coordinate(slice%params.slices_per_process,ptcl), new_coordinates[slice%params.slices_per_process], dist, params.box_size);
                            double r = sqrt(std::inner_product(dist.begin(), dist.end(), dist.begin(), 0.0));
                            new_action_local += potential_value_coulomb(dist, kVec, paths.charge[ptcl], paths.charge[ptcl1], params);
#pragma omp critical
                            new_distances.push_back(std::tuple<std::pair<int,int>,std::vector<double>,double>(std::pair<int,int>(paths.get_coordinate_key(slice%params.slices_per_process, ptcl),paths.get_coordinate_key(slice%params.slices_per_process, ptcl1)),  dist, r));
                        }
                    }
                    if(i != 0 && i != params.Mbar){
                        dist = std::vector<double>(params.dimensions, 0);
                        new_action_local += potential_value_coulomb(dist, kVec, paths.charge[ptcl1], paths.charge[ptcl1], params);
                    }
                    old_action_p[slice%params.slices_per_process] += old_action_local;
                    new_action_p[slice%params.slices_per_process] += new_action_local;
                    break;
            }
        }
    }
    if(id != 0){
        MPI_Send(&old_action_p[0], params.slices_per_process, MPI_DOUBLE, 0,0,local_comm);
        MPI_Send(&new_action_p[0], params.slices_per_process, MPI_DOUBLE, 0,1,local_comm);
    }
    else{
        ac_re = false;
        std::vector<double> new_action;
        std::vector<double> old_action;
        new_action.reserve(params.total_slices);
        old_action.reserve(params.total_slices);
        old_action.insert(old_action.end(),old_action_p.begin(),old_action_p.end());
        new_action.insert(new_action.end(),new_action_p.begin(),new_action_p.end());
        for(int j = 1; j < params.num_workers; ++j){
            MPI_Recv(&old_action_p[0], params.slices_per_process, MPI_DOUBLE,j,0,local_comm,MPI_STATUS_IGNORE);
            old_action.insert(old_action.end(),old_action_p.begin(),old_action_p.end());
            MPI_Recv(&new_action_p[0], params.slices_per_process, MPI_DOUBLE,j,1,local_comm,MPI_STATUS_IGNORE);
            new_action.insert(new_action.end(),new_action_p.begin(),new_action_p.end());
        }
        double na = 0;
        double oa = 0;
        for(int j = 0; j < params.Mbar; ++j){
            int slice = (start_slice+j)%params.total_slices;
            na += params.tau/2*(new_action[slice]+new_action[(slice+1)%params.total_slices]);
            oa += params.tau/2*(old_action[slice]+old_action[(slice+1)%params.total_slices]);
        }
        if(rng.randnormed(1) < (sigma_I/sigma_Z)*exp(-(na-oa)))
            ac_re = true;
    }
}

int Swap_Head::attempt(int &id, Parameters &params, Paths &paths, RNG &rng){
    if(!params.worm_on) return 0;
    Moves::attempt(id, params, paths, rng);
    int start_slice = params.worm_head.first;
    int swap_slice = (start_slice - params.Mbar + params.total_slices) % params.total_slices;
    keep_going = true;
    sigma_I = 0;
    sigma_Z = 0;
    std::vector<double> start(params.dimensions,0);
    std::vector<double> end(params.dimensions,0);
    choice = 0;
    if(start_slice >= params.my_start && start_slice <= params.my_end){
        int sslice = start_slice%params.slices_per_process;
        start = paths.get_coordinate(sslice, params.worm_head.second);
        if(swap_slice >= params.my_start && swap_slice <= params.my_end){
            int eslice = swap_slice%params.slices_per_process;
            std::vector<int> possible_swaps = paths.get_nearest_neighbors(eslice, start);
            if(!possible_swaps.empty()){
                for(int poss = possible_swaps.size() - 1; poss >= 0; --poss){
                    if(paths.charge[possible_swaps[poss]] != paths.charge[params.worm_tail.second])
                        possible_swaps.erase(possible_swaps.begin() + poss);
                }
            }
            if(possible_swaps.size() == 0)
                keep_going = false;
            else{
                std::vector<double> rho0s(possible_swaps.size(),0);
                std::vector<double> dist(params.dimensions,0);
                for(auto &i : possible_swaps){
                    distance(start, paths.get_coordinate(eslice, i), dist, params.box_size);
                    double r2 = inner_product(dist.begin(), dist.end(), dist.begin(),0.0);
                    rho0s[&i - &possible_swaps[0]] = pow((4* M_PI*params.lambda*params.Mbar*params.tau),-params.dimensions/2.)*exp(-1.0/(4.0*params.lambda*params.Mbar*params.tau)*r2);
                }
                for (auto &r : rho0s)
                    sigma_I += r;
                double total = 0;
                double rn = rng.randnormed(sigma_I);
                for(auto &r: rho0s){
                    total += r;
                    if(rn <= total){
                        choice = possible_swaps[&r - &rho0s[0]];
                        break;
                    }
                }
                if(paths.broken[choice]){
                    int for_part = choice;
                    if(start_slice - params.Mbar < 0)
                        for_part = paths.forward_connects[choice];
                    if(for_part == -1 || (for_part == params.worm_tail.second && start_slice > params.worm_tail.first))
                        keep_going = false;
                }
                if(keep_going){
                    int choice_back = choice;
                    if(start_slice - params.Mbar < 0)
                            choice_back = paths.forward_connects[choice];
                    std::vector<double> cb_loc(params.dimensions,0);
                    end = paths.get_coordinate(eslice, choice);
                    cb_loc = paths.get_coordinate(sslice,choice_back);
                    if((keep_going = paths.are_neighbors(cb_loc, end))){
                        std::vector<int> possible_swaps_2 = paths.get_nearest_neighbors(eslice, cb_loc);
                        if(!possible_swaps_2.empty()){
                            for(int poss = possible_swaps_2.size() - 1; poss >= 0; --poss){
                                if(paths.charge[possible_swaps_2[poss]] != paths.charge[params.worm_tail.second])
                                    possible_swaps_2.erase(possible_swaps_2.begin() + poss);
                            }
                        }
                        std::vector<double> rho0s_2(possible_swaps_2.size(),0);
                        for(auto &i : possible_swaps_2){
                            distance(cb_loc, paths.get_coordinate(eslice, i), dist, params.box_size);
                            double r2 = inner_product(dist.begin(),dist.end(), dist.begin(),0.0);
                            rho0s_2[&i - &possible_swaps_2[0]] = pow((4* M_PI*params.lambda*params.Mbar*params.tau),-params.dimensions/2.)*exp(-1.0/(4.0*params.lambda*params.Mbar*params.tau)*r2);
                        }
                        for (auto &r : rho0s_2)
                            sigma_Z += r;
                    }
                }
            }
        }
        else{
            MPI_Send(&start[0], params.dimensions, MPI_DOUBLE,swap_slice/params.slices_per_process, 0, local_comm);
            MPI_Recv(&keep_going, 1, MPI_INT, swap_slice/params.slices_per_process, 1, local_comm, MPI_STATUS_IGNORE);
            if(keep_going){
                MPI_Recv(&keep_going, 1, MPI_INT, swap_slice/params.slices_per_process, 2, local_comm, MPI_STATUS_IGNORE);
                if(keep_going){
                    MPI_Recv(&choice, 1, MPI_INT, swap_slice/params.slices_per_process, 3, local_comm, MPI_STATUS_IGNORE);
                    int choice_back = choice;
                    if(start_slice - params.Mbar <  0){
                        choice_back = paths.forward_connects[choice];
                    }
                    std::vector<double> cb_loc(params.dimensions,0);
                    cb_loc = paths.get_coordinate(sslice,choice_back);
                    MPI_Send(&cb_loc[0], params.dimensions, MPI_DOUBLE,swap_slice/params.slices_per_process, 4, local_comm);
                }
            }
        }
    }
    else if(swap_slice >= params.my_start && swap_slice <= params.my_end){
        std::vector<double> tail(params.dimensions,0);
        MPI_Recv(&tail[0], params.dimensions, MPI_DOUBLE, start_slice/params.slices_per_process, 0, local_comm, MPI_STATUS_IGNORE);
        int slice = swap_slice%params.slices_per_process;
        std::vector<int> possible_swaps = paths.get_nearest_neighbors(slice, tail);
        if(!possible_swaps.empty()){
            for(int poss = possible_swaps.size() - 1; poss >= 0; --poss){
                if(paths.charge[possible_swaps[poss]] != paths.charge[params.worm_tail.second])
                    possible_swaps.erase(possible_swaps.begin() + poss);
            }
        }
        if(possible_swaps.size() == 0)
            keep_going = false;
        MPI_Send(&keep_going, 1, MPI_INT, start_slice/params.slices_per_process, 1, local_comm);
        if(keep_going){
            std::vector<double> rho0s(possible_swaps.size(),0);
            std::vector<double> dist(params.dimensions,0);
            for(auto &i : possible_swaps){
                distance(tail, paths.get_coordinate(slice, i), dist, params.box_size);
                double r2 = inner_product(dist.begin(),dist.end(), dist.begin(),0.0);
                rho0s[&i - &possible_swaps[0]] = pow((4* M_PI*params.lambda*params.Mbar*params.tau),-params.dimensions/2.)*exp(-1.0/(4.0*params.lambda*params.Mbar*params.tau)*r2);
            }
            for (auto &r : rho0s)
                sigma_I += r;
            double total = 0;
            double rn = rng.randnormed(sigma_I);
            for(auto &r : rho0s){
                total += r;
                if(rn <= total){
                    choice = possible_swaps[&r - &rho0s[0]];
                    break;
                }
            }
            if(paths.broken[choice]){
                int for_part = choice;
                if(start_slice - params.Mbar < 0)
                    for_part = paths.forward_connects[choice];
                if(for_part == -1 || (for_part == params.worm_tail.second && start_slice > params.worm_tail.first))
                    keep_going = false;
            }
            MPI_Send(&keep_going, 1, MPI_INT, start_slice/params.slices_per_process, 2, local_comm);
            if(keep_going){
                MPI_Send(&choice, 1, MPI_INT, start_slice/params.slices_per_process, 3, local_comm);
                end = paths.get_coordinate(slice, choice);
                std::vector<double> cb(params.dimensions,0);
                MPI_Recv(&cb[0], params.dimensions, MPI_DOUBLE, start_slice/params.slices_per_process, 4, local_comm, MPI_STATUS_IGNORE);
                if((keep_going = paths.are_neighbors(cb, end))){
                    std::vector<int> possible_swaps_2 = paths.get_nearest_neighbors(slice, cb);
                    if(!possible_swaps_2.empty()){
                        for(int poss = possible_swaps_2.size() - 1; poss >= 0; --poss){
                            if(paths.charge[possible_swaps_2[poss]] != paths.charge[params.worm_tail.second])
                                possible_swaps_2.erase(possible_swaps_2.begin() + poss);
                        }
                    }
                    std::vector<double> rho0s_2(possible_swaps_2.size(),0);
                    for(auto &i : possible_swaps_2){
                        distance(end, paths.get_coordinate(slice, i), dist, params.box_size);
                        double r2 = inner_product(dist.begin(),dist.end(), dist.begin(),0.0);
                        rho0s_2[&i - &possible_swaps_2[0]] = pow((4* M_PI*params.lambda*params.Mbar*params.tau),-params.dimensions/2.)*exp(-1.0/(4.0*params.lambda*params.Mbar*params.tau)*r2);
                    }
                    for (auto &r : rho0s_2)
                        sigma_Z += r;
                }
            }
        }
    }
    MPI_Bcast(&keep_going, 1, MPI_INT, swap_slice/params.slices_per_process, local_comm);
    if(!keep_going) return 1;
    MPI_Bcast(&choice, 1, MPI_INT, swap_slice/params.slices_per_process, local_comm);
    MPI_Bcast(&sigma_I, 1, MPI_DOUBLE, swap_slice/params.slices_per_process, local_comm);
    MPI_Bcast(&sigma_Z, 1, MPI_DOUBLE, swap_slice/params.slices_per_process, local_comm);
    MPI_Bcast(&end[0], params.dimensions, MPI_DOUBLE, swap_slice/params.slices_per_process, local_comm);
    new_coordinates.resize(params.slices_per_process);
    new_distances.reserve(params.Mbar*params.particles);
    if(swap_slice >= params.my_start && swap_slice <= params.my_end)
        new_coordinates[swap_slice%params.slices_per_process] = end;
    for(int i = 1; i < params.Mbar; ++i){
        int slice = (start_slice-i + params.total_slices)%params.total_slices;
        int slicep1 = (start_slice-i+1+params.total_slices)%params.total_slices;
        if(slice >= params.my_start && slice <= params.my_end){
            if(slicep1 >= params.my_start && slicep1 <= params.my_end){
                if(i == 1){
                    new_coordinates[slicep1%params.slices_per_process] = start;
                }
                start = new_coordinates[slicep1%params.slices_per_process];
            }
            else
                MPI_Recv(&start[0], params.dimensions, MPI_DOUBLE, slicep1/params.slices_per_process, i, local_comm, MPI_STATUS_IGNORE);
            average_loc_weighted(start, end, new_coordinates[slice%params.slices_per_process], params.box_size, 1, params.Mbar-i);
            double tau_w = params.tau*params.Mbar-i;
            double width = sqrt(2.*params.lambda/(1./params.tau+1./tau_w));
            for(auto &j : new_coordinates[slice%params.slices_per_process])
                j += rng.randgaussian(width);
            put_in_box(new_coordinates[slice%params.slices_per_process], params.box_size);
        }
        else{
            if(slicep1 >= params.my_start && slicep1 <= params.my_end){
                if(i == 1)
                    new_coordinates[slicep1%params.slices_per_process] = start;
                MPI_Send(&new_coordinates[slicep1%params.slices_per_process][0], params.dimensions,MPI_DOUBLE,slice/params.slices_per_process, i, local_comm);
            }
        }
    }
    check(id, params, paths, rng);
    MPI_Bcast(&ac_re, 1, MPI_INT,0, local_comm);
    ++num_attempts;
    if(ac_re){
        ++num_accepts;
        paths.swap_into_head(params, choice, params.Mbar, new_coordinates);
        paths.update_separations(new_distances);
    }
    new_coordinates.clear();
    new_distances.clear();
    return 0;
}

void Swap_Head::check(int &id, Parameters &params, Paths &paths, RNG &rng){
    std::vector<double> old_action_p(params.slices_per_process,0);
    std::vector<double> new_action_p(params.slices_per_process,0);
    std::vector<double> dist(params.dimensions,0);
    std::vector<int> kVec(params.dimensions,0);
    if(params.potential == 2){
        params.alpha = sqrt(M_PI)*pow(params.particles/params.volume2,1/6.);
        params.coulcut = sqrt(params.p)/params.alpha;
        params.coulcut2 = pow(params.coulcut,2);
        params.kcut = 2.*params.p/params.coulcut;
        params.kcut2 = pow(params.kcut,2);
        params.nmax = floor(params.coulcut/params.box_size);
        params.kmax = ceil(params.kcut/(2.*M_PI/params.box_size));
    }
    int start_slice = params.worm_head.first;
    for(int i = 0; i <= params.Mbar; ++i){
        int ptcl1 = choice;
        int slice = (start_slice-i+params.total_slices)%params.total_slices;
        if(slice >= params.my_start && slice <= params.my_end){
            if(start_slice-params.Mbar < 0 && start_slice-i >= 0)
                ptcl1 = paths.forward_connects[ptcl1];
            double old_action_local = 0;
            double new_action_local = 0;
            switch(params.potential){
                case 0:
                {
                    new_action_p[slice%params.slices_per_process] = potential_value_harmonic(inner_product(new_coordinates[slice%params.slices_per_process].begin(), new_coordinates[slice%params.slices_per_process].end(), new_coordinates[slice%params.slices_per_process].begin(),0.0));
                    old_action_p[slice%params.slices_per_process] = potential_value_harmonic(inner_product(paths.get_coordinate(slice%params.slices_per_process, ptcl1).begin(), paths.get_coordinate(slice%params.slices_per_process, ptcl1).end(), paths.get_coordinate(slice%params.slices_per_process, ptcl1).begin(), 0.0));
                    break;
                }
                case 1:
#pragma omp parallel for firstprivate(dist) reduction(+:new_action_local,old_action_local)
                    for(int ptcl = 0; ptcl < params.particles; ++ptcl){
                        if(!(ptcl == params.worm_head.second && slice < params.worm_head.first) && !(ptcl == params.worm_tail.second && slice > params.worm_tail.first)){
                            if(i != 0 && ptcl != ptcl1){
                                old_action_local += potential_value_aziz(paths.get_separation(slice%params.slices_per_process, ptcl, ptcl1));
                                distance(paths.get_coordinate(slice%params.slices_per_process,ptcl), new_coordinates[slice%params.slices_per_process], dist, params.box_size);
                                double r = sqrt(std::inner_product(dist.begin(), dist.end(), dist.begin(), 0.0));
                                new_action_local += potential_value_aziz(r);
#pragma omp critical
                                new_distances.push_back(std::tuple<std::pair<int,int>,std::vector<double>,double>(std::pair<int,int>(paths.get_coordinate_key(slice%params.slices_per_process, ptcl),paths.get_coordinate_key(slice%params.slices_per_process, ptcl1)),  dist, r));
                            }
                            else if(i == 0){
                                if(ptcl != ptcl1)
                                    old_action_local += potential_value_aziz(paths.get_separation(slice%params.slices_per_process, ptcl, ptcl1));
                                distance(paths.get_coordinate(slice%params.slices_per_process,ptcl), new_coordinates[slice%params.slices_per_process], dist, params.box_size);
                                if(ptcl != params.worm_head.second)
                                    new_action_local += potential_value_aziz(sqrt(std::inner_product(dist.begin(), dist.end(), dist.begin(), 0.0)));
                            }
                        }
                    }
                    old_action_p[slice%params.slices_per_process] += old_action_local;
                    new_action_p[slice%params.slices_per_process] += new_action_local;
                    break;
                case 2:
#pragma omp parallel for firstprivate(dist) reduction(+:new_action_local,old_action_local)
                    for(int ptcl = 0; ptcl < params.particles; ++ptcl){
                        if(!(ptcl == params.worm_head.second && slice < params.worm_head.first) && !(ptcl == params.worm_tail.second && slice > params.worm_tail.first)){
                            old_action_local += potential_value_coulomb(paths.get_separation_vector(slice%params.slices_per_process, ptcl, ptcl1), kVec, paths.charge[ptcl], paths.charge[ptcl1], params);
                            distance(paths.get_coordinate(slice%params.slices_per_process,ptcl), new_coordinates[slice%params.slices_per_process], dist, params.box_size);
                            double r = sqrt(std::inner_product(dist.begin(), dist.end(), dist.begin(), 0.0));
                            new_action_local += potential_value_coulomb(dist, kVec, paths.charge[ptcl], paths.charge[ptcl1], params);
#pragma omp critical
                            new_distances.push_back(std::tuple<std::pair<int,int>,std::vector<double>,double>(std::pair<int,int>(paths.get_coordinate_key(slice%params.slices_per_process, ptcl),paths.get_coordinate_key(slice%params.slices_per_process, ptcl1)),  dist, r));
                        }
                    }
                    if(i != 0 && i != params.Mbar){
                        dist = std::vector<double>(params.dimensions, 0);
                        new_action_local += potential_value_coulomb(dist, kVec, paths.charge[ptcl1], paths.charge[ptcl1], params);
                    }
                    old_action_p[slice%params.slices_per_process] += old_action_local;
                    new_action_p[slice%params.slices_per_process] += new_action_local;
                    break;
            }
        }
    }
    if(id != 0){
        MPI_Send(&old_action_p[0], params.slices_per_process, MPI_DOUBLE, 0,0,local_comm);
        MPI_Send(&new_action_p[0], params.slices_per_process, MPI_DOUBLE, 0,1,local_comm);
    }
    else{
        ac_re = false;
        std::vector<double> new_action;
        std::vector<double> old_action;
        new_action.reserve(params.total_slices);
        old_action.reserve(params.total_slices);
        old_action.insert(old_action.end(),old_action_p.begin(),old_action_p.end());
        new_action.insert(new_action.end(),new_action_p.begin(),new_action_p.end());
        for(int j = 1; j < params.num_workers; ++j){
            MPI_Recv(&old_action_p[0], params.slices_per_process, MPI_DOUBLE,j,0,local_comm,MPI_STATUS_IGNORE);
            old_action.insert(old_action.end(),old_action_p.begin(),old_action_p.end());
            MPI_Recv(&new_action_p[0], params.slices_per_process, MPI_DOUBLE,j,1,local_comm,MPI_STATUS_IGNORE);
            new_action.insert(new_action.end(),new_action_p.begin(),new_action_p.end());
        }
        double na = 0;
        double oa = 0;
        for(int j = 0; j < params.Mbar; ++j){
            int slice = (start_slice-j+params.total_slices)%params.total_slices;
            na += params.tau/2*(new_action[slice]+new_action[(slice+1)%params.total_slices]);
            oa += params.tau/2*(old_action[slice]+old_action[(slice+1)%params.total_slices]);
        }
        if(rng.randnormed(1) < (sigma_I/sigma_Z)*exp(-(na-oa)))
            ac_re = true;
    }
}
