//
//  lookuptable.hpp
//  BEC-monopoles
//
//  Created by Adith Ramamurti on 11/16/15.
//  Copyright © 2015 Adith Ramamurti. All rights reserved.
//

#ifndef lookuptable_hpp
#define lookuptable_hpp

#include "utility.h"

class Separation_Table{
    
public:
    
    Separation_Table(double boxsize = -1){
        util = boost::shared_ptr<Utility>(new Utility(0));
        this->boxsize = boxsize;
    };
    ~Separation_Table(){};
    
    void add_bead(int slice, size_t key, std::vector<double> location){
        if(slice >= beads_map.size())
            beads_map.resize(slice+1);
        beads_map[slice].insert(std::pair<size_t, std::vector<double> >(key, location));
        // std::cout << key << "\t"<< slice << std::endl;
        
        update_bead(slice,key,location);
    }
    
    void remove_bead(int slice, size_t key){
        auto bead = beads_map[slice].find(key);
        for(boost::unordered_map<size_t, std::vector<double> >::iterator it = beads_map[slice].begin(); it != beads_map[slice].end(); it++){
            if(it != bead){
                std::pair<size_t,size_t> bead_pair_id_1 = std::pair<size_t,size_t>(bead->first,it->first);
                std::pair<size_t,size_t> bead_pair_id_2 = std::pair<size_t,size_t>(it->first,bead->first);
                bead_pair_sep.erase(bead_pair_id_1);
                bead_pair_sep.erase(bead_pair_id_2);
            }
        }
        
        // std::cout << key << std::endl;
        beads_map[slice].erase(key);
    }
    
    
    void update_bead(int slice, size_t key, dVector new_loc){
        auto bead = beads_map[slice].find(key);
        bead->second = new_loc;

        for(boost::unordered_map<size_t, std::vector<double> >::iterator it = beads_map[slice].begin(); it != beads_map[slice].end(); it++){
            if(it != bead){
                ddVector bead_pair;
                std::pair<size_t, size_t> bead_pair_id_1 = std::pair<size_t,size_t>(bead->first,it->first);
                std::pair<size_t, size_t> bead_pair_id_2= std::pair<size_t,size_t>(it->first,bead->first);
                bead_pair.push_back(bead->second);
                bead_pair.push_back(it->second);
                dVector dist = util->dist(bead_pair,boxsize);
                
                auto it2 = bead_pair_sep.find(bead_pair_id_1);
                if(it2 == bead_pair_sep.end()){
                    bead_pair_sep.insert(std::pair<std::pair<size_t,size_t>, dVector>(bead_pair_id_1,dist));
                }
                else
                   bead_pair_sep.find(bead_pair_id_1)->second = dist;

                dVector neg_dist;
                std::transform(dist.begin(), dist.end(), std::back_inserter(neg_dist), std::negate<double>());
                   
                it2 = bead_pair_sep.find(bead_pair_id_2);
                if(it2 == bead_pair_sep.end()){
                    bead_pair_sep.insert(std::pair<std::pair<size_t, size_t>, dVector>(bead_pair_id_2,neg_dist));
                }
                else
                    bead_pair_sep.find(bead_pair_id_2)->second = neg_dist;
            }
        }
    }
    
    dVector get_separation(std::pair<size_t, size_t> bead_pair){
        auto it = bead_pair_sep.find(bead_pair);
        return it->second;
    }
    
    void set_old_maps(){
        old_beads_map = beads_map;
        old_bead_pair_sep = bead_pair_sep;
    }
    void revert_maps(){
        beads_map = old_beads_map;
        bead_pair_sep = old_bead_pair_sep;
    }
    
    
private:
    
    std::vector<boost::unordered_map<size_t, std::vector<double> > > beads_map;
    boost::unordered_map<std::pair<size_t,size_t>, std::vector<double> > bead_pair_sep;
    
    std::vector<boost::unordered_map<size_t, std::vector<double> > > old_beads_map;
    boost::unordered_map<std::pair<size_t,size_t>, std::vector<double> > old_bead_pair_sep;
    
    
    boost::shared_ptr<Utility> util;
    
    double boxsize;
    
};

template<typename A, typename B, typename C>
struct boost::hash<boost::tuple<A,B,C> >{
    size_t operator()(const boost::tuple<A,B,C> &t) const{
        size_t seed = 0;
        boost::hash_combine(seed, t.template get<0>());
        boost::hash_combine(seed, t.template get<1>());
        boost::hash_combine(seed, t.template get<2>());
        return seed;
    }
};

template<class T>

class Neighbor_Table{
    
public:
    
    Neighbor_Table(double boxsize, int ndim){
        this->box_size = boxsize;
        num_grid = int(ceil(box_size/4));
        grid_step = boxsize/num_grid;
        this->ndim = ndim;
        util = boost::shared_ptr<Utility>(new Utility(0));
        
    };
    
    ~Neighbor_Table(){};
    
    void set_up_neighbor_table(){
        size_t grid_key;
        for(int i = 0; i < num_grid; i++){
            std::vector<int> grid(ndim);
            grid[0]=i;
            if(ndim>1)
                for(int j = 0; j < num_grid; j++){
                    grid[1]=j;
                    if(ndim>2)
                        for(int k = 0; k < num_grid; k++){
                            grid[2]=k;
                            grid_key = triple_hash(boost::tuple<int,int,int>(grid[0],grid[1],grid[2]));
                            iiVector neighboring_boxes = get_shifted_box(grid);
                            std::vector<size_t> neighbor_box_keys(0);
                            for(iiVector::iterator it = neighboring_boxes.begin(); it != neighboring_boxes.end(); it++){
                                neighbor_box_keys.push_back(triple_hash(boost::tuple<int,int,int>((*it)[0],(*it)[1],(*it)[2])));
                            }
                            grid_neigbors.insert(std::pair<size_t, std::vector<size_t> >(grid_key, neighbor_box_keys));
                            
                        }
                    else{
                        grid_key = pair_hash(std::pair<int, int>(grid[0],grid[1]));
                        iiVector neighboring_boxes = get_shifted_box(grid);
                        std::vector<size_t> neighbor_box_keys(0);
                        for(iiVector::iterator it = neighboring_boxes.begin(); it != neighboring_boxes.end(); it++){
                            neighbor_box_keys.push_back(pair_hash(std::pair<int,int>((*it)[0],(*it)[1])));
                        }
                        grid_neigbors.insert(std::pair<size_t, std::vector<size_t> >(grid_key, neighbor_box_keys));
                    }
                }
            else{
                iiVector neighboring_boxes = get_shifted_box(grid);
                std::vector<size_t> neighbor_box_keys(0);
                for(iiVector::iterator it = neighboring_boxes.begin(); it != neighboring_boxes.end(); it++){
                    neighbor_box_keys.push_back((*it)[0]);
                }
                grid_neigbors.insert(std::pair<size_t, std::vector<size_t> >(grid[0], neighbor_box_keys));
            }
        }
    }
    
    iiVector get_shifted_box(std::vector<int> grid_box){
        iiVector shifts;
        
        switch(ndim){
            case 1:
                for(int i = -1; i <= 1; i++)
                    if(i != 0){
                        std::vector<int> shift;
                        shift.push_back(i);
                        shifts.push_back(shift);
                    }
                
                break;
            case 2:
                for(int i = -1; i <= 1; i++)
                    for(int j = -1; j <= 1; j++)
                        if(i != 0 || j != 0){
                            std::vector<int> shift;
                            shift.push_back(i);
                            shift.push_back(j);
                            shifts.push_back(shift);
                        }
                break;
            case 3:
                for(int i = -1; i <= 1; i++)
                    for(int j = -1; j <= 1; j++)
                        for(int k = -1; k <= 1; k++)
                            if(i != 0 || j != 0 || k != 0){
                                std::vector<int> shift;
                                shift.push_back(i);
                                shift.push_back(j);
                                shift.push_back(k);
                                shifts.push_back(shift);
                            }
                break;
        }
        
        
        iiVector neighboring_boxes;
        
        for(iiVector::iterator it = shifts.begin(); it != shifts.end(); it++){
            iVector shift = *it;
            iVector result;
            result.reserve(ndim);
            std::transform(grid_box.begin(), grid_box.end(), shift.begin(), std::back_inserter(result), std::plus<int>());
            for(iVector::iterator it2 = result.begin(); it2 != result.end(); it2++)
                if(*it2 < 0)
                    *it2 = num_grid-1;
            neighboring_boxes.push_back(result);
        }
        
        return neighboring_boxes;
    }
    
    
    void add_bead(T reference){
        dVector data = reference->data;
        int col = reference->column_number;

        if(data.size()!=0){
            iVector grid_num;
            for(int i = 0; i < ndim; i++){
                grid_num.push_back(util->per_bound_cond((int)(data[i]/grid_step), num_grid));
            }
            size_t grid_key;
            switch(ndim){
                case 1:
                    grid_key = grid_num[0];
                    break;
                case 2:
                    grid_key = pair_hash(std::pair<int, int>(grid_num[0],grid_num[1]));
                    break;
                case 3:
                    grid_key = triple_hash(boost::tuple<int,int,int>(grid_num[0],grid_num[1],grid_num[2]));
                    break;
            }
            if(grid_beads.size() <= col)
                grid_beads.resize(col+1);
            auto it = grid_beads[col].find(grid_key);
            if(it == grid_beads[col].end()){
                std::vector<size_t> keys;
                keys.push_back(reference->key);
                grid_beads[col].insert(std::pair<size_t, std::vector<size_t> >(grid_key, keys));
            }
            else{
                it->second.push_back(reference->key);
            }
        }
    }
    
    
    void update_bead(T reference){
        dVector data = reference->data;
        dVector old_data = reference->old_data;
       
        iVector old_grid_num(0);

        if(old_data.size()!=0){
            for(int i = 0; i < ndim; i++)
                old_grid_num.push_back(util->per_bound_cond((int)(old_data[i]/grid_step),num_grid));
        }
        
        iVector new_grid_num(0);
        for(int i = 0; i < ndim; i++)
            new_grid_num.push_back(util->per_bound_cond((int)(data[i]/grid_step),num_grid));
        
        if(std::equal(old_grid_num.begin(), old_grid_num.end(), new_grid_num.begin()) && old_grid_num.size() != 0){
            return;
        }
        
        int col = reference->column_number;
        
        size_t new_grid_key = 0;
        switch(ndim){
            case 1:
                new_grid_key = new_grid_num[0];
                break;
            case 2:
                new_grid_key = pair_hash(std::pair<int, int>(new_grid_num[0],new_grid_num[1]));
                break;
            case 3:
                new_grid_key = triple_hash(boost::tuple<int,int,int>(new_grid_num[0],new_grid_num[1],new_grid_num[2]));
                break;
        }
        
        size_t old_grid_key = 0;
        
        if(old_data.size()!=0){
            switch(ndim){
                case 1:
                    old_grid_key = old_grid_num[0];
                    break;
                case 2:
                    old_grid_key = pair_hash(std::pair<int, int>(old_grid_num[0],old_grid_num[1]));
                    break;
                case 3:
                    old_grid_key = triple_hash(boost::tuple<int,int,int>(old_grid_num[0],old_grid_num[1],old_grid_num[2]));
                    break;
            }
            
            auto it = grid_beads[col].find(old_grid_key);
            if(it != grid_beads[col].end()){
                auto it2 = std::find(it->second.begin(), it->second.end(), reference->key);
                it->second.erase(it2);
            }
        }
        
        auto it2 = grid_beads[col].find(new_grid_key);
        if(it2 == grid_beads[col].end()){
            std::vector<size_t> keys;
            keys.push_back(reference->key);
            grid_beads[col].insert(std::pair<size_t, std::vector<size_t> >(new_grid_key, keys));
        }
        else{
            it2->second.push_back(reference->key);
        }
        
    }
    
    
    void remove_bead(T reference){
        dVector data = reference->data;
        iVector grid_num;
        
        size_t grid_key = 0;
        int col = reference->column_number;

        if(data.size() != 0){
            for(int i = 0; i < ndim; i++)
                grid_num.push_back(util->per_bound_cond((int)(data[i]/grid_step),num_grid));
            
            switch(ndim){
                case 1:
                    grid_key = grid_num[0];
                    break;
                case 2:
                    grid_key = pair_hash(std::pair<int, int>(grid_num[0],grid_num[1]));
                    break;
                case 3:
                    grid_key = triple_hash(boost::tuple<int,int,int>(grid_num[0],grid_num[1],grid_num[2]));
                    break;
            }
        
 
        auto it = grid_beads[col].find(grid_key);
        if(it != grid_beads[col].end()){
            auto it2 = std::find(it->second.begin(), it->second.end(), reference->key);
            it->second.erase(it2);
        }
        }
    }
    
    bool check_bead(T reference){
        dVector data = reference->data;
        if(data.size() != 0){
            iVector grid_num;
            for(int i = 0; i < ndim; i++)
                grid_num.push_back(util->per_bound_cond((int)(data[i]/grid_step),num_grid));
            
            size_t grid_key;
            int col = reference->column_number;
            switch(ndim){
                case 1:
                    grid_key = grid_num[0];
                    break;
                case 2:
                    grid_key = pair_hash(std::pair<int, int>(grid_num[0],grid_num[1]));
                    break;
                case 3:
                    grid_key = triple_hash(boost::tuple<int,int,int>(grid_num[0],grid_num[1],grid_num[2]));
                    break;
            }
            
            auto it = grid_beads[col].find(grid_key);
            std::vector<size_t> keys = it->second;
            auto it2 = std::find(keys.begin(), keys.end(), reference->key);
            if(it2 != keys.end())
                return true;
            else
                return false;
        }
        return true;
    }
    
    size_t get_bead_grid_num(T reference){
        dVector data = reference->data;
        if(data.size() != 0){
            iVector grid_num;
            for(int i = 0; i < ndim; i++)
                grid_num.push_back(util->per_bound_cond((int)(data[i]/grid_step),num_grid));
            
            size_t grid_key = 0;
            switch(ndim){
                case 1:
                    grid_key = grid_num[0];
                    break;
                case 2:
                    grid_key = pair_hash(std::pair<int, int>(grid_num[0],grid_num[1]));
                    break;
                case 3:
                    grid_key = triple_hash(boost::tuple<int,int,int>(grid_num[0],grid_num[1],grid_num[2]));
                    break;
            }
            return grid_key;
        }
        return 0;
    }
    
    void set_old_table(){
        old_grid_beads = grid_beads;
    }
    
    void reset_old_table(){
        grid_beads = old_grid_beads;
    }
    
    std::vector<size_t> get_neighboring_beads(T reference, int col){
        size_t grid_key = get_bead_grid_num(reference);
        std::vector<size_t> all_grid_boxes = grid_neigbors.find(grid_key)->second;
        all_grid_boxes.push_back(grid_key);
        std::vector<size_t> bead_keys(0);
        for(std::vector<size_t>::iterator it = all_grid_boxes.begin(); it != all_grid_boxes.end(); it++){
            std::vector<size_t> addtl_bead_keys = get_beads_in_box(*it,col);
            bead_keys.insert(bead_keys.end(),addtl_bead_keys.begin(),addtl_bead_keys.end());
        }
        return bead_keys;
    }
    
    std::vector<size_t> get_beads_in_box(size_t grid_key, int col){
        boost::unordered_map<size_t, std::vector<size_t> > grid_bead_col = grid_beads[col];
        auto it =grid_bead_col.find(grid_key);
        std::vector<size_t> bead_keys(0);
        if(it != grid_bead_col.end())
            bead_keys = it->second;
        return bead_keys;
    }
    
    
private:
    
    int num_grid;
    double box_size;
    double grid_step;
    int ndim;
    
    std::vector<boost::unordered_map<size_t, std::vector<size_t> > > grid_beads;
    std::vector<boost::unordered_map<size_t, std::vector<size_t> > > old_grid_beads;
    
    boost::unordered_map<size_t, std::vector<size_t> > grid_neigbors;
    
    boost::hash<boost::tuple<int, int, int> > triple_hash;
    boost::hash<std::pair<int, int> > pair_hash;
    boost::shared_ptr<Utility> util;
    
    
};

#endif /* lookuptable_hpp */
