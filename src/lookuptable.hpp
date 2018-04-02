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

//hashing methods for two and three integers: ((n1+n2)*(n1+n2+1))/2+n2 for two ints, and recursively for three
inline size_t pair_hash(std::pair<int, int> pair){
    return ((pair.first+pair.second)*(pair.first+pair.second+1))/2+pair.second;
}

inline size_t triple_hash(std::tuple<int, int, int> tuple){
    return pair_hash(std::make_pair(pair_hash(std::make_pair(std::get<0>(tuple), std::get<1>(tuple))), std::get<2>(tuple)));
}

/***This class holds a neighbor table for the particles. It divides n-D space into a grid, and stores which particles are in neighboring grid boxes for quick lookup.***/
class Neighbor_Table{
    
public:
    Neighbor_Table(Parameters &params){ //sets box size and neighbor table grid size from parameters
        if(params.box_size != -1){
            box_size = params.box_size;
            box_start = 0;
        }
        else{
            box_size = 30;
            box_start = -box_size/2;
        }
        if(params.grid_size != 0){
            num_grid = ceil(box_size/params.grid_size);
            num_grid = std::max(num_grid,3);
        }
        else
            num_grid = 3;
        grid_step = box_size/num_grid;
        cell_members.resize(params.slices_per_process);
        set_up_neighbor_table(params);
    }
    
    ~Neighbor_Table(){}
    
    //Sets up the neighbor table. Each grid box is assighned a unique key, and the neighboring grid boxes for any grid box are stored in a vector
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
                            for(std::vector<std::vector<int> >::iterator it = neighboring_boxes.begin(); it != neighboring_boxes.end(); ++it)
                                neighbor_box_keys.push_back(triple_hash(std::make_tuple((*it)[0],(*it)[1],(*it)[2])));
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
    
    //helper method to find neighboring boxes
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
            for(std::vector<int>::iterator it2 = result.begin(); it2 != result.end(); ++it2){
                if(*it2 < 0)
                    *it2 = num_grid-1;
                if(*it2 >= num_grid)
                    *it2 = *it2-num_grid;
            }
            neighboring_boxes.push_back(result);
        }
        return neighboring_boxes;
    }
    
    //adds a bead to the table (takes bead slice, bead key, and bead position)
    void add_bead(int slice, int key, const std::vector<double>& position){
        size_t grid_key = get_grid_key(position);
        auto found = cell_members[slice].find(grid_key);
        if(found != cell_members[slice].end())
            (found->second).push_back(key);
        else
            cell_members[slice].insert({grid_key,{key}});
    }
    
    //removes a bead from the table
    void remove_bead(int slice, int key, const std::vector<double>& position){
        if(position.size() != 0){
            size_t grid_key = get_grid_key(position);
            auto found = cell_members[slice].find(grid_key);
            if(found != cell_members[slice].end())
                (found->second).erase(std::remove(found->second.begin(), found->second.end(), key), found->second.end());
        }
    }
    
    //updates a bead's position in the grid
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
    
    //gets the nearest neighbors of a bead given slice and position. returns neighboring keys.
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
    
    //gets beads in the same box as the given bead
    std::vector<int> get_cell_mates(int slice, const std::vector<double>& position){
        size_t grid_key = get_grid_key(position);
        std::vector<int> neighboring_particles(0);
        auto found = cell_members[slice].find(grid_key);
        if(found != cell_members[slice].end())
            return found->second;
        return neighboring_particles;
    }
    
    //checks whether two beads are nearest grid neighbors
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
    
    //gets the key of the grid box given a position vector
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
    
    
    /******* Checker methods ************/
    
    //prints the whole map to the console (grid key: {bead keys})
    void print_map(int slice){
        for(auto &it1 : cell_members[slice]){
            std::cout << it1.first << ":\t";
            for(auto &it2 : it1.second)
                std::cout << it2 << " ";
            std::cout << std::endl;
        }
    }
    
    //gets the number of beads in the map at any slice
    int get_members_size(int slice){
        int size = 0;
        for(auto &i : cell_members[slice]){
            size += i.second.size();
        }
        return size;
    }
    
    //checks whether a bead is stored in the correct box
    bool is_in_box(int slice, int key, std::vector<double>& pos){
        size_t grid_key = get_grid_key(pos);
        auto cell = cell_members[slice].find(grid_key);
        if(cell != cell_members[slice].end()){
            auto vit = std::find(cell->second.begin(), cell->second.end(), key);
            if(vit != cell->second.end())
                return true;
        }
        return false;
    }
    
private:
    int num_grid;
    double box_size;
    double box_start;
    double grid_step;
    std::unordered_map<size_t, std::vector<size_t> > grid_neighbors;
    std::vector<std::unordered_map<size_t, std::vector<int> > > cell_members;
};

/*** This class stores the separations between all beads in a particular slice ***/
class Separation_Table{
    
public:
    Separation_Table(Parameters &params){
        locations.resize(params.slices_per_process);
        box_size = params.box_size;
    }
    
    ~Separation_Table(){}
    
    //adds a bead to the table
    void add_bead(int slice, int key, const std::vector<double>& position, bool update = true){
        locations[slice].insert({key,position});
        if(update)
            calculate_separations(slice, key);
    }
    
    //removes a bead from the table
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
                potentials.erase(key_fwd);
                potentials.erase(key_bwd);
            }
            locations[slice].erase(i);
        }
    }
    
    //updates the location of a bead
    void update_location(int slice, int key, const std::vector<double>& new_position, bool update = true){
        auto found = locations[slice].find(key);
        if(found != locations[slice].end())
            found->second = new_position;
        if(update)
            calculate_separations(slice, key);
    }
    
    //updates separations between beads given a set of pairs of keys, new distance vectors, and new distances
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
    
    //updates the potential values between beads given pair of keys and the potential value
    void update_potentials(std::vector<std::tuple<std::pair<int,int>, double> >& new_potentials){
        for(auto &n : new_potentials){
            size_t key_fwd = pair_hash(std::get<0>(n));
            size_t key_bwd = pair_hash(std::pair<int, int>(std::get<0>(n).second,std::get<0>(n).first));
            potentials[key_fwd] = std::get<1>(n);
            potentials[key_bwd] = std::get<1>(n);
        }
    }
    
    //calculates the separations between one bead and all other beads
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
    
    /*** Getter methods ***/
    
    std::vector<double>& get_distance_vector(int key1, int key2){
        size_t key = pair_hash(std::pair<int, int>(key1,key2));
        return distance_vectors[key];
    }
    
    double& get_distance(int key1, int key2){
        size_t key = pair_hash(std::pair<int, int>(key1,key2));
        return distances[key];
    }
    
    double& get_potential(int key1, int key2){
        size_t key = pair_hash(std::pair<int, int>(key1,key2));
        return potentials[key];
    }

    
private:
    std::vector<std::unordered_map<int, std::vector<double> > > locations;
    std::unordered_map<size_t, std::vector<double> > distance_vectors;
    std::unordered_map<size_t, double> distances;
    std::unordered_map<size_t, double> potentials;
    double box_size;
};


/*** This class stores separations between beads in the imaginary time direction (across slices) ***/
class Kinetic_Separation_Table{
    
public:
    Kinetic_Separation_Table(Parameters &params){
        start_locations.resize(params.slices_per_process);
        end_locations.resize(params.slices_per_process);
        box_size = params.box_size;
    }
    
    ~Kinetic_Separation_Table(){}
    
    void add_bead_start(int slice, int key, const std::vector<double>& position){
        start_locations[slice].insert({key,position});
    }
    
    void add_bead_end(int slice, int key, const std::vector<double>& position){
        end_locations[slice].insert({key,position});
    }

    //removes a bead from an end slice
    void remove_bead_end(int slice, int key){
        auto i = end_locations[slice].find(key);
        if(i != end_locations[slice].end()){
            for(auto& j : start_locations[slice]){
                size_t key_fwd = pair_hash(std::pair<int, int>(i->first,j.first));
                size_t key_bwd = pair_hash(std::pair<int, int>(j.first,i->first));
                distances.erase(key_fwd);
                distances.erase(key_bwd);
            }
            end_locations[slice].erase(i);
        }
    }
    
    //removes a bead from a start slice
    void remove_bead_start(int slice, int key){
        auto i = start_locations[slice].find(key);
        if(i != start_locations[slice].end()){
            for(auto& j : end_locations[slice]){
                size_t key_fwd = pair_hash(std::pair<int, int>(i->first,j.first));
                size_t key_bwd = pair_hash(std::pair<int, int>(j.first,i->first));
                distances.erase(key_fwd);
                distances.erase(key_bwd);
            }
            start_locations[slice].erase(i);
        }
    }
    
    //updates the location of a start bead
    void update_location_start(int slice, int key, const std::vector<double>& new_position, bool update = true){
        auto found = start_locations[slice].find(key);
        if(found != start_locations[slice].end())
            found->second = new_position;
        if(update)
            calculate_separations_start(slice, key);
    }
    
    //updates the location of an end bead
    void update_location_end(int slice, int key, const std::vector<double>& new_position, bool update = true){
        auto found = end_locations[slice].find(key);
        if(found != end_locations[slice].end())
            found->second = new_position;
        if(update)
            calculate_separations_end(slice, key);
    }

    
    // NOTE THAT THE SEPARATIONS STORED ARE SQUARE SEPARATIONS
    
    //calculates distance from the start particle to all end particles
    void calculate_separations_start(int slice, int ptcl){
        auto i = start_locations[slice].find(ptcl);
        if(i != start_locations[slice].end()){
            for(auto& j : end_locations[slice]){
                size_t key_fwd = pair_hash(std::pair<int, int>(i->first,j.first));
                size_t key_bwd = pair_hash(std::pair<int, int>(j.first,i->first));
                std::vector<double> dist(i->second.size());
                distance(i->second, j.second, dist, box_size);
                distances[key_fwd] = inner_product(dist.begin(), dist.end(), dist.begin(), 0.0);
                distances[key_bwd] = distances[key_fwd];
            }
        }
    }
    
    //calcules distances from the end particle to all beginning
    void calculate_separations_end(int slice, int ptcl){
        auto i = end_locations[slice].find(ptcl);
        if(i != end_locations[slice].end()){
            for(auto& j : start_locations[slice]){
                size_t key_fwd = pair_hash(std::pair<int, int>(i->first,j.first));
                size_t key_bwd = pair_hash(std::pair<int, int>(j.first,i->first));
                std::vector<double> dist(i->second.size());
                distance(i->second, j.second, dist, box_size);
                distances[key_fwd] = inner_product(dist.begin(), dist.end(), dist.begin(), 0.0);
                distances[key_bwd] = distances[key_fwd];
            }
        }
    }
    
    //gets the distance between any two beads
    double get_distance(int key1, int key2){
        size_t key = pair_hash(std::pair<int, int>(key1,key2));
        return distances[key];
    }
    
private:
    std::vector<std::unordered_map<int, std::vector<double> > > start_locations;
    std::vector<std::unordered_map<int, std::vector<double> > > end_locations;
    std::unordered_map<size_t, double> distances;
    double box_size;
};



#endif /* lookuptable_hpp */

