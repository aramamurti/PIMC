//
//  paths.hpp
//
//
//  Created by Adith Ramamurti on 9/12/16.
//
//

#ifndef paths_hpp
#define paths_hpp

#include <stdio.h>
#include <mpi.h>
#include <iomanip>

#include "utility.hpp"
#include "parameters.hpp"
#include "lookuptable.hpp"

inline int positive_modulo(int i, int n) {
    return (n + (i % n)) % n;
}

class Paths{
private:
    std::vector<std::vector<std::vector<double> > > coordinate_slices;
    std::vector<std::vector<int> > coordinate_keys;
    std::unordered_map<int, std::pair<int, int> > key_finder;
    Neighbor_Table nt;
    Separation_Table st;
    RNG rng;
    int bead_counter;
    MPI_Comm local_comm;
    
public:
    std::vector<int> forward_connects;
    std::vector<int> backward_connects;
    std::vector<bool> broken;
    std::vector<int> charge;
    int broken_worldlines;
    
    Paths(int &id, Parameters &params, MPI_Comm &local) : nt(params), st(params){
        local_comm = local;
        rng.seed(id);
        set_up(id, params);
    }
    
    void set_up(int &id, Parameters &params){
        std::vector<double> locations(params.particles*params.dimensions);
        if(id == 0){
            for(auto &i : locations){
                if(params.box_size != -1)
                    i = rng.randnormed(params.box_size);
                else
                    i = 1-rng.randnormed(2);
            }
        }
        MPI_Bcast(&locations[0], (int)locations.size(), MPI_DOUBLE, 0, local_comm);
        forward_connects.resize(params.particles);
        backward_connects.resize(params.particles);
        charge.resize(params.particles,0);
        broken.resize(params.particles, false);
        broken_worldlines = 0;
        for(int i = 0; i < params.particles; ++i){
            forward_connects[i] = i;
            backward_connects[i] = i;
        }
        if(params.charged)
            for(int i = 0; i < params.particles; ++i)
                charge[i] = (i%params.charges)*2-1;
        bead_counter = 1;
        coordinate_slices.resize(params.slices_per_process);
        coordinate_keys.resize(params.slices_per_process);
        int slice_counter = 0;
        for(auto &slice: coordinate_slices){
            slice.resize(params.particles);
            coordinate_keys[&slice - &coordinate_slices[0]].resize(params.particles);
            int counter = 0;
            for(auto &particle : slice){
                particle.resize(params.dimensions);
                for(auto &dim : particle){
                    dim = locations[counter];
                    ++counter;
                }
                coordinate_keys[&slice - &coordinate_slices[0]][&particle - &slice[0]] = bead_counter;
                key_finder.insert({bead_counter, std::pair<int, int>(&slice - &coordinate_slices[0],&particle - &slice[0])});
                nt.add_bead(slice_counter, bead_counter, particle);
                st.add_bead(slice_counter, bead_counter, particle);
                ++bead_counter;
            }
            ++slice_counter;
        }
        params.worm_on = false;
        params.worm_head.first = -1;
        params.worm_head.second = -1;
        params.worm_tail.first = -1;
        params.worm_tail.second = -1;
    }
    
    /*******************************************************************************************************************************************
     
     WORM METHODS
     
     *******************************************************************************************************************************************/

    void worm_advance_tail(Parameters& params, const std::vector<double>& location, const std::vector<std::tuple<int, std::vector<double>, double> >& distances, int worm_start = 0, int chg = -10){
        if(params.worm_on == false){
            params.worm_on = true;
            params.worm_head.second = params.particles;
            params.worm_tail.second = params.particles;
            params.worm_head.first = worm_start;
            params.worm_tail.first = worm_start;
            ++params.particles;
            forward_connects.push_back(-1);
            backward_connects.push_back(-1);
            broken.push_back(true);
            if(chg != -10)
                charge.push_back(chg);
            else
                charge.push_back(rng.randint(params.charges)*2-1);
        }
        else if((++params.worm_tail.first) == params.total_slices){
            forward_connects[params.worm_tail.second] = params.particles;
            forward_connects.push_back(-1);
            backward_connects.push_back(params.worm_tail.second);
            charge.push_back(charge[params.worm_tail.second]);
            params.worm_tail.first = 0;
            params.worm_tail.second = params.particles;
            broken.push_back(true);
            ++params.particles;
        }
        if(params.worm_tail.first >= params.my_start && params.worm_tail.first <= params.my_end){
            int my_slice = params.worm_tail.first%params.slices_per_process;
            if(params.worm_tail.second >= coordinate_slices[my_slice].size()){
                coordinate_slices[my_slice].resize(params.worm_tail.second+1);
                coordinate_keys[my_slice].resize(params.worm_tail.second+1);
            }
            coordinate_slices[my_slice][params.worm_tail.second] = location;
            coordinate_keys[my_slice][params.worm_tail.second] = bead_counter;
            key_finder.insert({bead_counter, std::pair<int, int>(my_slice,params.worm_tail.second)});
            nt.add_bead(my_slice, bead_counter, location);
            st.add_bead(my_slice, bead_counter, location, false);
            std::vector<std::tuple<std::pair<int, int>, std::vector<double>, double> > seps;
            seps.reserve(distances.size()+1);
            seps.push_back(std::tuple<std::pair<int, int>, std::vector<double>, double>(std::pair<int,int>(bead_counter,bead_counter), std::vector<double>(params.dimensions,0), 0));
            for(auto &sep : distances)
                seps.push_back(std::tuple<std::pair<int, int>, std::vector<double>, double>(std::pair<int,int>(std::get<0>(sep), bead_counter), std::get<1>(sep), std::get<2>(sep)));
            update_separations(seps);
            ++bead_counter;
        }
        ++params.worm_length;
    }
    
    void worm_advance_head(Parameters& params, const std::vector<double>& location, const std::vector<std::tuple<int, std::vector<double>, double> >& distances, int worm_start = 0, int chg = -10){
        if(!params.worm_on){
            params.worm_on = true;
            params.worm_head.second = params.particles;
            params.worm_tail.second = params.particles;
            params.worm_head.first = worm_start;
            params.worm_tail.first = worm_start;
            ++params.particles;
            forward_connects.push_back(-1);
            backward_connects.push_back(-1);
            if(chg != -10)
                charge.push_back(chg);
            else
                charge.push_back(rng.randint(params.charges)*2-1);
            broken.push_back(true);
            ++broken_worldlines;
        }
        else if((--params.worm_head.first) == -1){
            backward_connects[params.worm_head.second] = params.particles;
            backward_connects.push_back(-1);
            forward_connects.push_back(params.worm_head.second);
            charge.push_back(charge[params.worm_head.second]);
            params.worm_head.first = params.total_slices - 1;
            params.worm_head.second = params.particles;
            broken.push_back(true);
            ++broken_worldlines;
            ++params.particles;
        }
        if(params.worm_head.first >= params.my_start && params.worm_head.first <= params.my_end){
            int my_slice = params.worm_head.first%params.slices_per_process;
            if(params.worm_head.second >= coordinate_slices[my_slice].size()){
                coordinate_slices[my_slice].resize(params.worm_head.second+1);
                coordinate_keys[my_slice].resize(params.worm_head.second+1);
            }
            coordinate_slices[my_slice][params.worm_head.second] = location;
            coordinate_keys[my_slice][params.worm_head.second] = bead_counter;
            key_finder.insert({bead_counter, std::pair<int, int>(my_slice,params.worm_head.second)});
            nt.add_bead(my_slice, bead_counter, location);
            st.add_bead(my_slice, bead_counter, location, false);
            std::vector<std::tuple<std::pair<int, int>, std::vector<double>, double> > seps;
            seps.reserve(distances.size()+1);
            seps.push_back(std::tuple<std::pair<int, int>, std::vector<double>, double>(std::pair<int,int>(bead_counter,bead_counter), std::vector<double>(params.dimensions,0), 0));
            for(auto &sep : distances)
                seps.push_back(std::tuple<std::pair<int, int>, std::vector<double>, double>(std::pair<int,int>(std::get<0>(sep),bead_counter), std::get<1>(sep), std::get<2>(sep)));
            update_separations(seps);
            ++bead_counter;
        }
        ++params.worm_length;
    }
    
    void worm_recede_tail(Parameters& params){
        if(!params.worm_on) return;
        if(params.worm_tail.first >= params.my_start && params.worm_tail.first <= params.my_end){
            int my_slice = params.worm_tail.first%params.slices_per_process;
            nt.remove_bead(my_slice, coordinate_keys[my_slice][params.worm_tail.second], coordinate_slices[my_slice][params.worm_tail.second]);
            st.remove_bead(my_slice, coordinate_keys[my_slice][params.worm_tail.second]);
            key_finder.erase(coordinate_keys[my_slice][params.worm_tail.second]);
            coordinate_keys[my_slice][params.worm_tail.second] = 0;
            coordinate_slices[my_slice][params.worm_tail.second].clear();
        }
        if(--params.worm_tail.first == -1 && params.worm_length > 1){
            params.worm_tail.first = params.total_slices - 1;
            int col_to_remove = params.worm_tail.second;
            params.worm_tail.second = backward_connects[params.worm_tail.second];
            forward_connects[params.worm_tail.second] = -1;
            if(params.worm_head.second > col_to_remove) --params.worm_head.second;
            if(params.worm_tail.second > col_to_remove) --params.worm_tail.second;
            for(int i = 0; i < backward_connects.size(); ++i){
                if(backward_connects[i] > col_to_remove)
                    --backward_connects[i];
                if(forward_connects[i] > col_to_remove)
                    --forward_connects[i];
            }
            backward_connects.erase(backward_connects.begin() + col_to_remove);
            forward_connects.erase(forward_connects.begin() + col_to_remove);
            broken.erase(broken.begin() + col_to_remove);
            charge.erase(charge.begin() + col_to_remove);
            --broken_worldlines;
            for(int slice = 0; slice < params.slices_per_process; ++slice){
                if(coordinate_slices[slice].size() > col_to_remove){
                    coordinate_slices[slice].erase(coordinate_slices[slice].begin() + col_to_remove);
                    coordinate_keys[slice].erase(coordinate_keys[slice].begin() + col_to_remove);
                }
                for(int col = 0; col < coordinate_keys[slice].size(); ++col){
                    if(key_finder[coordinate_keys[slice][col]].second > col_to_remove)
                        --key_finder[coordinate_keys[slice][col]].second;
                }
            }
            --params.particles;
        }
        if(--params.worm_length == 0){
            params.worm_on = false;
            int col_to_remove = params.worm_tail.second;
            for(int i = 0; i < backward_connects.size(); ++i){
                if(backward_connects[i] > col_to_remove)
                    --backward_connects[i];
                if(forward_connects[i] > col_to_remove)
                    --forward_connects[i];
            }
            backward_connects.erase(backward_connects.begin() + col_to_remove);
            forward_connects.erase(forward_connects.begin() + col_to_remove);
            broken.erase(broken.begin() + col_to_remove);
            charge.erase(charge.begin() + col_to_remove);
            --broken_worldlines;
            for(int slice = 0; slice < params.slices_per_process; ++slice){
                if(coordinate_slices[slice].size() > col_to_remove){
                    coordinate_slices[slice].erase(coordinate_slices[slice].begin() + col_to_remove);
                    coordinate_keys[slice].erase(coordinate_keys[slice].begin() + col_to_remove);
                }
                for(int col = 0; col < coordinate_keys[slice].size(); ++col){
                    if(key_finder[coordinate_keys[slice][col]].second > col_to_remove)
                        --key_finder[coordinate_keys[slice][col]].second;
                }
            }
            --params.particles;
            params.worm_head.first = -1;
            params.worm_head.second = -1;
            params.worm_tail.first = -1;
            params.worm_tail.second = -1;
        }
    }
    
    void worm_recede_head(Parameters& params){
        if(!params.worm_on) return;
        if(params.worm_head.first >= params.my_start && params.worm_head.first <= params.my_end){
            int my_slice = params.worm_head.first%params.slices_per_process;
            nt.remove_bead(my_slice, coordinate_keys[my_slice][params.worm_head.second], coordinate_slices[my_slice][params.worm_head.second]);
            st.remove_bead(my_slice, coordinate_keys[my_slice][params.worm_head.second]);
            key_finder.erase(coordinate_keys[my_slice][params.worm_head.second]);
            coordinate_keys[my_slice][params.worm_head.second] = 0;
            coordinate_slices[my_slice][params.worm_head.second].clear();
        }
        if(++params.worm_head.first == params.total_slices && params.worm_length > 1){
            params.worm_head.first = 0;
            int col_to_remove = params.worm_head.second;
            params.worm_head.second = forward_connects[params.worm_head.second];
            backward_connects[params.worm_head.second] = -1;
            if(params.worm_head.second > col_to_remove) --params.worm_head.second;
            if(params.worm_tail.second > col_to_remove) --params.worm_tail.second;
            for(int i = 0; i < backward_connects.size(); ++i){
                if(backward_connects[i] > col_to_remove)
                    --backward_connects[i];
                if(forward_connects[i] > col_to_remove)
                    --forward_connects[i];
            }
            backward_connects.erase(backward_connects.begin() + col_to_remove);
            forward_connects.erase(forward_connects.begin() + col_to_remove);
            broken.erase(broken.begin() + col_to_remove);
            charge.erase(charge.begin() + col_to_remove);
            --broken_worldlines;
            for(int slice = 0; slice < params.slices_per_process; ++slice){
                if(coordinate_slices[slice].size() > col_to_remove){
                    coordinate_slices[slice].erase(coordinate_slices[slice].begin() + col_to_remove);
                    coordinate_keys[slice].erase(coordinate_keys[slice].begin() + col_to_remove);
                }
                for(int col = 0; col < coordinate_keys[slice].size(); ++col){
                    if(key_finder[coordinate_keys[slice][col]].second > col_to_remove)
                        --key_finder[coordinate_keys[slice][col]].second;
                }
            }
            --params.particles;
        }
        if(--params.worm_length == 0){
            params.worm_on = false;
            int col_to_remove = params.worm_head.second;
            for(int i = 0; i < backward_connects.size(); ++i){
                if(backward_connects[i] > col_to_remove)
                    --backward_connects[i];
                if(forward_connects[i] > col_to_remove)
                    --forward_connects[i];
            }
            backward_connects.erase(backward_connects.begin() + col_to_remove);
            forward_connects.erase(forward_connects.begin() + col_to_remove);
            broken.erase(broken.begin() + col_to_remove);
            charge.erase(charge.begin() + col_to_remove);
            --broken_worldlines;
            for(int slice = 0; slice < params.slices_per_process; ++slice){
                if(coordinate_slices[slice].size() > col_to_remove){
                    coordinate_slices[slice].erase(coordinate_slices[slice].begin() + col_to_remove);
                    coordinate_keys[slice].erase(coordinate_keys[slice].begin() + col_to_remove);
                }
                for(int col = 0; col < coordinate_keys[slice].size(); ++col){
                    if(key_finder[coordinate_keys[slice][col]].second > col_to_remove)
                        --key_finder[coordinate_keys[slice][col]].second;
                }
            }
            --params.particles;
            params.worm_head.first = -1;
            params.worm_head.second = -1;
            params.worm_tail.first = -1;
            params.worm_tail.second = -1;
        }
    }
    
    void open_path(Parameters& params, int column, int tail_slice, int distance){
        if(params.worm_on) return;
        params.worm_on = true;
        params.worm_length = 0;
        int current_column = column;
        do{
            params.worm_length += params.total_slices;
            broken[current_column] = true;
            ++broken_worldlines;
            current_column = backward_connects[current_column];
        }while(current_column != column);
        if(tail_slice != params.total_slices - 1){
            backward_connects.push_back(backward_connects[column]);
            forward_connects[backward_connects.back()] = params.particles;
            broken.push_back(true);
            charge.push_back(charge[column]);
            ++broken_worldlines;
            backward_connects[column] = -1;
            forward_connects.push_back(-1);
            for(int slice = 0; slice <= tail_slice; ++slice){
                if(slice >= params.my_start && slice <= params.my_end){
                    int my_slice = slice%params.slices_per_process;
                    coordinate_slices[my_slice].push_back(coordinate_slices[my_slice][column]);
                    coordinate_keys[my_slice].push_back(coordinate_keys[my_slice][column]);
                    key_finder[coordinate_keys[my_slice][column]] = std::pair<int, int>(my_slice, coordinate_slices[my_slice].size()-1);
                    coordinate_slices[my_slice][column].clear();
                    coordinate_keys[my_slice][column] = 0;
                }
            }
            params.worm_tail.first = tail_slice;
            params.worm_tail.second = params.particles;
            ++params.particles;
            params.worm_head.first = tail_slice+1;
            params.worm_head.second = column;
            if(params.worm_head.first == params.total_slices){
                params.worm_head.first = 0;
                params.worm_head.second = forward_connects[column];
            }
        }
        else{
            params.worm_tail.first = tail_slice;
            params.worm_tail.second = column;
            params.worm_head.first = 0;
            params.worm_head.second = forward_connects[column];
            backward_connects[forward_connects[column]] = -1;
            forward_connects[column] = -1;
        }
        for(int n = 0; n < distance; ++n)
            worm_recede_head(params);
    }
    
    void close_worm(Parameters& params, const std::vector<std::vector<double> >& new_coordinates, const std::vector<std::vector<std::tuple<int, std::vector<double>, double> > >& distances){
        if(!params.worm_on) return;
        while(params.worm_tail.first != positive_modulo(params.worm_head.first - 1, params.total_slices)){
            if(!new_coordinates.empty())
                worm_advance_tail(params, new_coordinates[(params.worm_tail.first+1)%params.slices_per_process], distances[(params.worm_tail.first+1)%params.slices_per_process]);
            else
                worm_advance_tail(params, std::vector<double>(params.dimensions, 0), std::vector<std::tuple<int, std::vector<double>, double> >());
        }
        int orig_worm_head = params.worm_head.first;
        while(params.worm_head.first != 0){
            if(params.worm_head.first >= params.my_start && params.worm_head.first <= params.my_end){
                int my_slice = params.worm_head.first%params.slices_per_process;
                if(coordinate_slices[my_slice].size() <= params.worm_tail.second){
                    coordinate_slices[my_slice].resize(params.worm_tail.second+1);
                    coordinate_keys[my_slice].resize(params.worm_tail.second+1);
                }
                coordinate_slices[my_slice][params.worm_tail.second] = coordinate_slices[my_slice][params.worm_head.second];
                coordinate_slices[my_slice][params.worm_head.second].clear();
                coordinate_keys[my_slice][params.worm_tail.second] = coordinate_keys[my_slice][params.worm_head.second];
                coordinate_keys[my_slice][params.worm_head.second] = 0;
                key_finder[coordinate_keys[my_slice][params.worm_tail.second]] = std::pair<int, int>(my_slice, params.worm_tail.second);
            }
            params.worm_tail.first = ++params.worm_tail.first%params.total_slices;
            params.worm_head.first = ++params.worm_head.first%params.total_slices;
        }
        if(orig_worm_head != 0){
            forward_connects[params.worm_tail.second] = forward_connects[params.worm_head.second];
            backward_connects[forward_connects[params.worm_tail.second]] = params.worm_tail.second;
            int col_to_remove = params.worm_head.second;
            if(params.worm_tail.second > col_to_remove) --params.worm_tail.second;
            for(int i = 0; i < backward_connects.size(); ++i){
                if(backward_connects[i] > col_to_remove)
                    --backward_connects[i];
                if(forward_connects[i] > col_to_remove)
                    --forward_connects[i];
            }
            backward_connects.erase(backward_connects.begin() + col_to_remove);
            forward_connects.erase(forward_connects.begin() + col_to_remove);
            broken.erase(broken.begin() + col_to_remove);
            charge.erase(charge.begin() + col_to_remove);
            for(int slice = 0; slice < params.slices_per_process; ++slice){
                if(coordinate_slices[slice].size() > col_to_remove){
                    coordinate_slices[slice].erase(coordinate_slices[slice].begin() + col_to_remove);
                    coordinate_keys[slice].erase(coordinate_keys[slice].begin() + col_to_remove);
                }
                for(int col = 0; col < coordinate_keys[slice].size(); ++col){
                    if(key_finder[coordinate_keys[slice][col]].second > col_to_remove)
                        --key_finder[coordinate_keys[slice][col]].second;
                }
            }
            --params.particles;
        }
        else{
            forward_connects[params.worm_tail.second] = params.worm_head.second;
            backward_connects[params.worm_head.second] = params.worm_tail.second;
        }
        broken_worldlines = 0;
        broken.clear();
        broken.resize(params.particles, false);
        params.worm_length = 0;
        params.worm_on = false;
        params.worm_head.first = -1;
        params.worm_head.second = -1;
        params.worm_tail.first = -1;
        params.worm_tail.second = -1;
    }
    
    void swap_into_head(Parameters& params, int column, int distance, const std::vector<std::vector<double> >& new_coordinates = std::vector<std::vector<double> >()){
        if(!params.worm_on) return;
        int new_head_col = column;
        if(params.worm_head.first-distance < 0){
            column = forward_connects[column];
            new_head_col = column;
        }
        int slice = (params.worm_head.first - 1);
        while(slice != -1){
            if(slice >= params.my_start && slice <= params.my_end){
                int my_slice = slice%params.slices_per_process;
                if(coordinate_slices[my_slice].size() <= params.worm_head.second){
                    coordinate_slices[my_slice].resize(params.worm_head.second+1);
                    coordinate_keys[my_slice].resize(params.worm_head.second+1);
                }
                coordinate_slices[my_slice][params.worm_head.second] = coordinate_slices[my_slice][column];
                coordinate_slices[my_slice][column].clear();
                coordinate_keys[my_slice][params.worm_head.second] = coordinate_keys[my_slice][column];
                coordinate_keys[my_slice][column] = 0;
                key_finder[coordinate_keys[my_slice][params.worm_head.second]] = std::pair<int, int>(my_slice, params.worm_head.second);
            }
            --slice;
        }
        if(params.worm_head.first-distance < 0){
            column = backward_connects[column];
            backward_connects[forward_connects[column]] = -1;
            backward_connects[params.worm_head.second] = column;
            forward_connects[column] = params.worm_head.second;
        }
        else{
            backward_connects[params.worm_head.second] = backward_connects[column];
            forward_connects[backward_connects[params.worm_head.second]] = params.worm_head.second;
            backward_connects[column] = -1;
        }
        if(!new_coordinates.empty()){
            slice = positive_modulo(params.worm_head.first - 1, params.total_slices);
            int col = params.worm_head.second;
            if(slice == params.total_slices - 1) col = backward_connects[col];
            for(int i = 0; i < distance-1; ++i){
                if(slice >= params.my_start && slice <= params.my_end){
                    int my_slice = slice%params.slices_per_process;
                    nt.update_grid(my_slice, coordinate_keys[my_slice][col], coordinate_slices[my_slice][col], new_coordinates[my_slice]);
                    st.update_location(my_slice, coordinate_keys[my_slice][col], new_coordinates[my_slice],false);
                    coordinate_slices[my_slice][col] = new_coordinates[my_slice];
                }
                --slice;
                if(slice == -1){
                    slice = params.total_slices - 1;
                    col = backward_connects[col];
                }
            }
        }
        params.worm_head.second = new_head_col;
        broken_worldlines = 0;
        broken.clear();
        broken.resize(params.particles, false);
        params.worm_length = 0;
        column = params.worm_head.second;
        while(column != -1){
            broken[column] = true;
            ++broken_worldlines;
            column = forward_connects[column];
            params.worm_length += params.total_slices;
        }
        params.worm_length -= (params.worm_head.first + (params.total_slices - params.worm_tail.first - 1));
    }
    
    void swap_into_tail(Parameters& params, int column, int distance, const std::vector<std::vector<double> >& new_coordinates = std::vector<std::vector<double> >()){
        if(!params.worm_on) return;
        int new_tail_col = column;
        if(params.worm_tail.first+distance >= params.total_slices){
            column = backward_connects[column];
            new_tail_col = column;
        }
        int slice = (params.worm_tail.first + 1);
        while(slice != params.total_slices){
            if(slice >= params.my_start && slice <= params.my_end){
                int my_slice = slice%params.slices_per_process;
                if(coordinate_slices[my_slice].size() <= params.worm_tail.second){
                    coordinate_slices[my_slice].resize(params.worm_tail.second+1);
                    coordinate_keys[my_slice].resize(params.worm_tail.second+1);
                }
                coordinate_slices[my_slice][params.worm_tail.second] = coordinate_slices[my_slice][column];
                coordinate_slices[my_slice][column].clear();
                coordinate_keys[my_slice][params.worm_tail.second] = coordinate_keys[my_slice][column];
                coordinate_keys[my_slice][column] = 0;
                key_finder[coordinate_keys[my_slice][params.worm_tail.second]] = std::pair<int, int>(my_slice, params.worm_tail.second);
            }
            ++slice;
        }
        if(params.worm_tail.first+distance >= params.total_slices){
            column = forward_connects[column];
            forward_connects[backward_connects[column]] = -1;
            forward_connects[params.worm_tail.second] = column;
            backward_connects[column] = params.worm_tail.second;
        }
        else{
            forward_connects[params.worm_tail.second] = forward_connects[column];
            backward_connects[forward_connects[params.worm_tail.second]] = params.worm_tail.second;
            forward_connects[column] = -1;
        }
        if(!new_coordinates.empty()){
            slice = (params.worm_tail.first + 1)%params.total_slices;
            int col = params.worm_tail.second;
            if(slice == 0) col = forward_connects[col];
            for(int i = 0; i < distance-1; ++i){
                if(slice >= params.my_start && slice <= params.my_end){
                    int my_slice = slice%params.slices_per_process;
                    nt.update_grid(my_slice, coordinate_keys[my_slice][col], coordinate_slices[my_slice][col], new_coordinates[my_slice]);
                    st.update_location(my_slice, coordinate_keys[my_slice][col], new_coordinates[my_slice], false);
                    coordinate_slices[my_slice][col] = new_coordinates[my_slice];
                }
                ++slice;
                if(slice == params.total_slices){
                    slice = 0;
                    col = forward_connects[col];
                }
            }
        }
        params.worm_tail.second = new_tail_col;
        broken_worldlines = 0;
        broken.clear();
        broken.resize(params.particles, false);
        params.worm_length = 0;
        column = params.worm_tail.second;
        while(column != -1){
            broken[column] = true;
            ++broken_worldlines;
            column = backward_connects[column];
            params.worm_length += params.total_slices;
        }
        params.worm_length -= (params.worm_head.first + (params.total_slices - params.worm_tail.first - 1));
    }
    
    /*******************************************************************************************************************************************
     
     GETTER & SETTER METHODS
     
     *******************************************************************************************************************************************/
    
    void set_coordinate(int slice, int column, const std::vector<double>& location, bool update = true){
        nt.update_grid(slice, coordinate_keys[slice][column], coordinate_slices[slice][column], location);
        st.update_location(slice, coordinate_keys[slice][column], location, update);
        coordinate_slices[slice][column] = location;
    }
    
    void update_separations(std::vector<std::tuple<std::pair<int,int>, std::vector<double>, double> >& new_distances){
        st.update_separations(new_distances);
    }
    
    std::vector<double>& get_coordinate(int slice, int column){
        return coordinate_slices[slice][column];
    }
    
    int get_coordinate_key(int slice, int column){
        return coordinate_keys[slice][column];
    }

    std::vector<int> get_nearest_neighbors(int slice, std::vector<double>& position){
        std::vector<int> keys = nt.get_nearest_neighbors(slice, position);
        std::vector<int> neighbors(keys.size());
        for(auto &key : keys){
            neighbors[&key - &keys[0]] = key_finder[key].second;
        }
        return neighbors;
    }
    
    bool are_neighbors(std::vector<double>& position1, std::vector<double>& position2){
        return nt.are_nearest_neighbors(position1, position2);
    }
    
    double get_separation(int slice, int ptcl1, int ptcl2){
        int key1 = coordinate_keys[slice][ptcl1];
        int key2 = coordinate_keys[slice][ptcl2];
        return st.get_distance(key1, key2);
    }
    
    std::vector<double>& get_separation_vector(int slice, int ptcl1, int ptcl2){
        int key1 = coordinate_keys[slice][ptcl1];
        int key2 = coordinate_keys[slice][ptcl2];
        return st.get_distance_vector(key1, key2);
    }
    
    void check_separations(Parameters& params){
        double s1 = 0;
        double s2 = 0;
        for(int s = 0; s < coordinate_slices.size(); ++s){
            for(int p1 = 0; p1 < coordinate_slices[s].size(); ++p1){
                for(int p2 = 0; p2 < coordinate_slices[s].size(); ++p2){
                    if(coordinate_slices[s][p1].size() != 0 && coordinate_slices[s][p2].size() != 0){
                        s1 += st.get_distance(coordinate_keys[s][p1], coordinate_keys[s][p2]);
                        std::vector<double> dist;
                        distance(coordinate_slices[s][p1], coordinate_slices[s][p2], dist, params.box_size);
                        s2 += sqrt(inner_product(dist.begin(),dist.end(),dist.begin(),0.0));
                    }
                }
            }
        }
        std::cout << s1 - s2 << std::endl;
        if(std::abs(s1 - s2) > 10E-10)
            getchar();
    }

    /*******************************************************************************************************************************************
     
     PRINT METHODS
     
     *******************************************************************************************************************************************/
    
    void print(int &id, Parameters &params){
        if(id == 0){
            int width = 26;
            std::cout << "FC:\t\t";
            for(auto &i : forward_connects)
                std::cout <<  std::setw(width) << std::left<< i ;
            std::cout << std::endl << "   \t\t";
            for(int i = 0; i < forward_connects.size(); ++i)
                std::cout << std::setw(width) << std::left<< "^";
            std::cout << std::endl << "   \t\t";
            for(int i = 0; i < forward_connects.size(); ++i)
                std::cout << std::setw(width)<< std::left<<"|";
            std::cout << std::endl;
            bool start = true;
            MPI_Send(&start, 1, MPI_INT, params.num_workers-1, 0, local_comm);
            MPI_Recv(&start, 1, MPI_INT, 1%params.num_workers,0,local_comm,MPI_STATUS_IGNORE);
            if(start){
                for(int i = params.slices_per_process-1; i >= 0; --i){
                    std::cout << i + params.my_start <<":\t";
                    for(auto &j : coordinate_slices[i]){
                        std::cout << "{";
                        for(int k = 0; k < j.size(); k++)
                            std::cout << std::fixed << std::setprecision(3) << j[k] <<" ";
                        if(j.size() == 0)
                            std::cout << std::setw(6*params.dimensions) <<" ";
                        std::cout << std::setw(6) << std::left <<  "} ";
                    }
                    std::cout << "||\n";
                }
                std::cout <<"   \t\t";
                for(int i = 0; i < backward_connects.size(); ++i)
                    std::cout << std::setw(width)<< std::left<<"|";
                std::cout << std::endl << "   \t\t";
                for(int i = 0; i < backward_connects.size(); ++i)
                    std::cout << std::setw(width) << std::left<< "v";
                std::cout << std::endl<< "BC:\t\t";
                for(auto &i : backward_connects)
                    std::cout <<  std::setw(width) << std::left<< i;
                std::cout << std::endl<< "BR:\t\t";
                for(int i = 0; i < broken.size(); i++)
                    std::cout <<  std::setw(width) << std::left<< broken[i];

                std::cout << std::endl;
                std::cout << std::endl;
            }
        }
        else{
            bool start = false;
            MPI_Recv(&start,1,MPI_INT, (id+1)%(params.num_workers),0,local_comm,MPI_STATUS_IGNORE);
            if(start)
                for(int i = params.slices_per_process-1; i >= 0; --i){
                    std::cout << i + params.my_start <<":\t";
                    for(auto &j : coordinate_slices[i]){
                        std::cout << "{";
                        for(int k = 0; k < j.size(); k++)
                            std::cout << std::fixed << std::setprecision(3) << j[k] <<" ";
                        if(j.size() == 0)
                            std::cout << std::setw(6*params.dimensions) <<" ";
                        std::cout << std::setw(6) << std::left <<  "} ";
                    }
                    std::cout << "||\n";
                }
            MPI_Send(&start, 1, MPI_INT, (id-1), 0, local_comm);
        }
    }
    
    void print_keys(int &id, Parameters &params){
        if(id == 0){
            int width = 16;
            std::cout << "FC:\t\t";
            for(auto &i : forward_connects)
                std::cout <<  std::setw(width) << std::left<< i ;
            std::cout << std::endl << "   \t\t";
            for(int i = 0; i < forward_connects.size(); ++i)
                std::cout << std::setw(width) << std::left<< "^";
            std::cout << std::endl << "   \t\t";
            for(int i = 0; i < forward_connects.size(); ++i)
                std::cout << std::setw(width)<< std::left<<"|";
            std::cout << std::endl;
            bool start = true;
            MPI_Send(&start, 1, MPI_INT, params.num_workers-1, 0, local_comm);
            MPI_Recv(&start, 1, MPI_INT, 1%params.num_workers,0,local_comm,MPI_STATUS_IGNORE);
            if(start){
                for(int i = params.slices_per_process-1; i >= 0; --i){
                    std::cout << i + params.my_start <<":\t";
                    for(auto &j : coordinate_keys[i]){
                        std::cout << "(" << std::setw(4) << j << ": "  <<  key_finder[j].first << ", " << key_finder[j].second << ")\t";
                    }
                    std::cout << "||\n";
                }
            }
            std::cout <<"   \t\t";
            for(int i = 0; i < backward_connects.size(); ++i)
                std::cout << std::setw(width)<< std::left<<"|";
            std::cout << std::endl << "   \t\t";
            for(int i = 0; i < backward_connects.size(); ++i)
                std::cout << std::setw(width) << std::left<< "v";
            std::cout << std::endl<< "BC:\t\t";
            for(auto &i : backward_connects)
                std::cout <<  std::setw(width) << std::left<< i;
            std::cout << std::endl<< "BR:\t\t";
            for(int i = 0; i < broken.size(); i++)
                std::cout <<  std::setw(width) << std::left<< broken[i];
            
            std::cout << std::endl;
            std::cout << std::endl;
        }
        else{
            bool start = false;
            MPI_Recv(&start,1,MPI_INT, (id+1)%(params.num_workers),0,local_comm,MPI_STATUS_IGNORE);
            if(start)
                for(int i = params.slices_per_process-1; i >= 0; --i){
                    std::cout << i + params.my_start <<":\t";
                    for(auto &j : coordinate_keys[i]){
                        std::cout << "(" << std::setw(4) << j << ": " <<  key_finder[j].first << ", " << key_finder[j].second << ")\t";
                    }
                    std::cout << "||\n";
                }
            MPI_Send(&start, 1, MPI_INT, (id-1), 0, local_comm);
        }
    }
    
    void print_nt(Parameters& params){
        for(int slice = 0; slice < params.total_slices; ++slice){
            std::cout << slice << std::endl;
            nt.print_map(slice);
            std::cout <<"\n" << std::endl;
        }
    }
};

#endif /* paths_hpp */
