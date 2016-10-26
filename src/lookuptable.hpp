//
//  lookuptable.hpp
//
//
//  Created by Adith Ramamurti on 9/19/16.
//
//

#ifndef lookuptable_hpp
#define lookuptable_hpp

#include "utility.hpp"
#include "parameters.hpp"

#include <unordered_map>


inline size_t pair_hash(std::pair<int, int> pair){
    return ((pair.first+pair.second)*(pair.first+pair.second+1))/2+pair.second;
}

inline size_t triple_hash(std::tuple<int, int, int> tuple){
    return pair_hash(std::make_pair(pair_hash(std::make_pair(std::get<0>(tuple), std::get<1>(tuple))), std::get<2>(tuple)));
}

class Neighbor_Table{
    
public:
    Neighbor_Table(Parameters &params){
        if(params.box_size != -1){
            box_size = params.box_size;
            box_start = 0;
        }
        else{
            box_size = 30;
            box_start = -box_size/2;
        }
        num_grid = 4;
        grid_step = box_size/num_grid;
        cell_members.resize(params.slices_per_process);
        set_up_neighbor_table(params);
    }
    
    ~Neighbor_Table(){}
    
    void set_up_neighbor_table(Parameters &params){
        size_t grid_key;
        for(int i = 0; i < num_grid; ++i){
            std::vector<int> grid(params.dimensions);
            grid[0]=i;
            if(params.dimensions>1)
                for(int j = 0; j < num_grid; ++j){
                    grid[1]=j;
                    if(params.dimensions>2)
                        for(int k = 0; k < num_grid; ++k){
                            grid[2]=k;
                            grid_key = triple_hash(std::make_tuple(grid[0],grid[1],grid[2]));
                            std::vector<std::vector<int> > neighboring_boxes = get_shifted_box(grid, params);
                            std::vector<size_t> neighbor_box_keys(0);
                            for(std::vector<std::vector<int> >::iterator it = neighboring_boxes.begin(); it != neighboring_boxes.end(); ++it){
                                neighbor_box_keys.push_back(triple_hash(std::make_tuple((*it)[0],(*it)[1],(*it)[2])));
                            }
                            std::sort(neighbor_box_keys.begin(), neighbor_box_keys.end());
                            auto last = std::unique(neighbor_box_keys.begin(), neighbor_box_keys.end());
                            neighbor_box_keys.erase(last, neighbor_box_keys.end());
                            grid_neighbors.insert(std::pair<size_t, std::vector<size_t> >(grid_key, neighbor_box_keys));
                        }
                    else{
                        grid_key = pair_hash(std::pair<int, int>(grid[0],grid[1]));
                        std::vector<std::vector<int> > neighboring_boxes = get_shifted_box(grid, params);
                        std::vector<size_t> neighbor_box_keys(0);
                        for(std::vector<std::vector<int> >::iterator it = neighboring_boxes.begin(); it != neighboring_boxes.end(); ++it){
                            neighbor_box_keys.push_back(pair_hash(std::pair<int,int>((*it)[0],(*it)[1])));
                        }
                        std::sort(neighbor_box_keys.begin(), neighbor_box_keys.end());
                        auto last = std::unique(neighbor_box_keys.begin(), neighbor_box_keys.end());
                        neighbor_box_keys.erase(last, neighbor_box_keys.end());
                        grid_neighbors.insert(std::pair<size_t, std::vector<size_t> >(grid_key, neighbor_box_keys));
                    }
                }
            else{
                std::vector<std::vector<int> > neighboring_boxes = get_shifted_box(grid, params);
                std::vector<size_t> neighbor_box_keys(0);
                for(std::vector<std::vector<int> >::iterator it = neighboring_boxes.begin(); it != neighboring_boxes.end(); ++it){
                    neighbor_box_keys.push_back((*it)[0]);
                }
                std::sort(neighbor_box_keys.begin(), neighbor_box_keys.end());
                auto last = std::unique(neighbor_box_keys.begin(), neighbor_box_keys.end());
                neighbor_box_keys.erase(last, neighbor_box_keys.end());
                grid_neighbors.insert(std::pair<size_t, std::vector<size_t> >(grid[0], neighbor_box_keys));
            }
        }
    }
    
    std::vector<std::vector<int> > get_shifted_box(std::vector<int>& grid_box, Parameters &params){
        std::vector<std::vector<int> > shifts;
        switch(params.dimensions){
            case 1:
                for(int i = -1; i <= 1; ++i)
                    if(i != 0){
                        std::vector<int> shift;
                        shift.push_back(i);
                        shifts.push_back(shift);
                    }
                break;
            case 2:
                for(int i = -1; i <= 1; ++i)
                    for(int j = -1; j <= 1; ++j)
                        if(i != 0 || j != 0){
                            std::vector<int> shift;
                            shift.push_back(i);
                            shift.push_back(j);
                            shifts.push_back(shift);
                        }
                break;
            case 3:
                for(int i = -1; i <= 1; ++i)
                    for(int j = -1; j <= 1; ++j)
                        for(int k = -1; k <= 1; ++k)
                            if(i != 0 || j != 0 || k != 0){
                                std::vector<int> shift;
                                shift.push_back(i);
                                shift.push_back(j);
                                shift.push_back(k);
                                shifts.push_back(shift);
                            }
                break;
        }
        std::vector<std::vector<int> > neighboring_boxes;
        for(std::vector<std::vector<int> >::iterator it = shifts.begin(); it != shifts.end(); ++it){
            std::vector<int> shift = *it;
            std::vector<int> result;
            result.reserve(params.dimensions);
            std::transform(grid_box.begin(), grid_box.end(), shift.begin(), std::back_inserter(result), std::plus<int>());
            for(std::vector<int>::iterator it2 = result.begin(); it2 != result.end(); ++it2)
                if(*it2 < 0)
                    *it2 = num_grid-1;
            neighboring_boxes.push_back(result);
        }
        return neighboring_boxes;
    }
    
    void add_bead(int slice, int key, const std::vector<double>& position){
        size_t grid_key = get_grid_key(position);
        auto found = cell_members[slice].find(grid_key);
        if(found != cell_members[slice].end())
            (found->second).push_back(key);
        else
            cell_members[slice].insert({grid_key,{key}});
    }
    
    void remove_bead(int slice, int key, const std::vector<double>& position){
        if(position.size() != 0){
            size_t grid_key = get_grid_key(position);
            auto found = cell_members[slice].find(grid_key);
            if(found != cell_members[slice].end())
                (found->second).erase(std::remove(found->second.begin(), found->second.end(), key), found->second.end());
        }
    }
    
    void update_grid(int slice, int key, const std::vector<double>& old_position, const std::vector<double>& new_position){
        size_t grid_key_old = get_grid_key(old_position);
        size_t grid_key_new = get_grid_key(new_position);
        if(grid_key_old == grid_key_new) return;
        auto found = cell_members[slice].find(grid_key_old);
        if(found != cell_members[slice].end())
            (found->second).erase(std::remove(found->second.begin(), found->second.end(), key), found->second.end());
        found = cell_members[slice].find(grid_key_new);
        if(found != cell_members[slice].end())
            (found->second).push_back(key);
        else
            cell_members[slice].insert({grid_key_new,{key}});
    }
    
    std::vector<int> get_nearest_neighbors(int slice, const std::vector<double>& position){
        size_t grid_key = get_grid_key(position);
        std::vector<size_t> neighbors = grid_neighbors.find(grid_key)->second;
        neighbors.push_back(grid_key);
        std::vector<int> neighboring_particles(0);
        for(auto &i : neighbors){
            auto found = cell_members[slice].find(i);
            if(found != cell_members[slice].end()){
                neighboring_particles.reserve(neighboring_particles.size() + found->second.size());
                neighboring_particles.insert(neighboring_particles.end(), found->second.begin(), found->second.end());
            }
        }
        return neighboring_particles;
    }
    
    bool are_nearest_neighbors(const std::vector<double>& position1, const std::vector<double>& position2){
        size_t grid_key_1 = get_grid_key(position1);
        size_t grid_key_2 = get_grid_key(position2);
        if(grid_key_1 == grid_key_2) return true;
        std::vector<size_t> neighbors = grid_neighbors.find(grid_key_1)->second;
        if(std::find(neighbors.begin(), neighbors.end(), grid_key_2) != neighbors.end())
            return true;
        else
            return false;
    }
    
    size_t get_grid_key(const std::vector<double>& position){
        std::vector<int> grid;
        for(auto &i : position)
            grid.push_back((i-box_start)/grid_step);
        if(grid.size() == 1){
            return grid[0];
        }
        else if(grid.size() == 2){
            return pair_hash(std::pair<int,int>(grid[0],grid[1]));
        }
        else{
            return triple_hash(std::make_tuple(grid[0],grid[1],grid[2]));
        }
    }
    
    void print_map(int slice){
        for(auto &it1 : cell_members[slice]){
            std::cout << it1.first << ":\t";
            for(auto &it2 : it1.second)
                std::cout << it2 << " ";
            std::cout << std::endl;
        }
    }
    
    int get_members_size(int slice){
        int size = 0;
        for(auto &i : cell_members[slice]){
            size += i.second.size();
        }
        return size;
    }
    
private:
    int num_grid;
    double box_size;
    double box_start;
    double grid_step;
    std::unordered_map<size_t, std::vector<size_t> > grid_neighbors;
    std::vector<std::unordered_map<size_t, std::vector<int> > > cell_members;
};

class Separation_Table{
    
public:
    Separation_Table(Parameters &params){
        locations.resize(params.slices_per_process);
        box_size = params.box_size;
    }
    
    ~Separation_Table(){}
    
    void add_bead(int slice, int key, const std::vector<double>& position, bool update = true){
        locations[slice].insert({key,position});
        if(update)
            calculate_separations(slice, key);
    }
    
    void remove_bead(int slice, int key){
        auto i = locations[slice].find(key);
        if(i != locations[slice].end()){
            for(auto& j : locations[slice]){
                size_t key_fwd = pair_hash(std::pair<int, int>(i->first,j.first));
                size_t key_bwd = pair_hash(std::pair<int, int>(j.first,i->first));
                distance_vectors.erase(key_fwd);
                distance_vectors.erase(key_bwd);
                distances.erase(key_fwd);
                distances.erase(key_bwd);
            }
            locations[slice].erase(i);
        }
    }
    
    void update_location(int slice, int key, const std::vector<double>& new_position, bool update = true){
        auto found = locations[slice].find(key);
        if(found != locations[slice].end())
            found->second = new_position;
        if(update)
            calculate_separations(slice, key);
    }
    
    void update_separations(std::vector<std::tuple<std::pair<int,int>, std::vector<double>, double> >& new_distances){
        for(auto &n : new_distances){
            size_t key_fwd = pair_hash(std::get<0>(n));
            size_t key_bwd = pair_hash(std::pair<int, int>(std::get<0>(n).second,std::get<0>(n).first));
            distance_vectors[key_fwd] = std::get<1>(n);
            distance_vectors[key_bwd].resize(0);
            std::transform(distance_vectors[key_fwd].begin(), distance_vectors[key_fwd].end(), std::back_inserter(distance_vectors[key_bwd]), std::negate<double>());
            distances[key_fwd] = std::get<2>(n);
            distances[key_bwd] = std::get<2>(n);
        }
    }
    
    void calculate_separations(int slice, int ptcl){
        auto i = locations[slice].find(ptcl);
        if(i != locations[slice].end()){
            for(auto& j : locations[slice]){
                size_t key_fwd = pair_hash(std::pair<int, int>(i->first,j.first));
                size_t key_bwd = pair_hash(std::pair<int, int>(j.first,i->first));
                std::vector<double> dist(i->second.size());
                distance(i->second, j.second, dist, box_size);
                distance_vectors[key_fwd] = dist;
                distance_vectors[key_bwd].resize(0);
                std::transform(distance_vectors[key_fwd].begin(), distance_vectors[key_fwd].end(), std::back_inserter(distance_vectors[key_bwd]), std::negate<double>());
                distances[key_fwd] = sqrt(inner_product(dist.begin(), dist.end(), dist.begin(), 0.0));
                distances[key_bwd] = distances[key_fwd];
            }
        }
    }
    
    std::vector<double>& get_distance_vector(int key1, int key2){
        size_t key = pair_hash(std::pair<int, int>(key1,key2));
        return distance_vectors[key];
    }
    
    double get_distance(int key1, int key2){
        size_t key = pair_hash(std::pair<int, int>(key1,key2));
        return distances[key];
    }
    
private:
    std::vector<std::unordered_map<int, std::vector<double> > > locations;
    std::unordered_map<size_t, std::vector<double> > distance_vectors;
    std::unordered_map<size_t, double> distances;
    double box_size;
};


#endif /* lookuptable_hpp */

