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
        update = true;
    };
    ~Separation_Table(){};
    
    void add_bead(int slice, size_t key, std::vector<double>& location){
        if(slice >= beads_map.size())
            beads_map.resize(slice+1);
        beads_map[slice].insert(std::pair<size_t, std::vector<double>& >(key, location));
        std::pair<size_t, size_t> bead_pair = std::pair<size_t,size_t>(key,key);
        bead_pair_sep[bead_pair] = std::vector<double>(location.size(),0);
    }
    
    void update_bead_separations(int slice, size_t key){
        boost::unordered_map<size_t, std::vector<double>& >::iterator bead;
        bead = beads_map[slice].find(key);
        
        for(boost::unordered_map<size_t, std::vector<double>& >::iterator it = beads_map[slice].begin(); it != beads_map[slice].end(); it++){
            if(it->first != bead->first){

                std::pair<size_t, size_t> bead_pair_id_1 = std::pair<size_t,size_t>(bead->first,it->first);
                std::pair<size_t, size_t> bead_pair_id_2= std::pair<size_t,size_t>(it->first,bead->first);

                if(update){
                    util->dist(bead->second, it->second, bead_pair_sep[bead_pair_id_1], boxsize);
                    bead_pair_sep[bead_pair_id_2].resize(0);
                    std::transform(bead_pair_sep[bead_pair_id_1].begin(), bead_pair_sep[bead_pair_id_1].end(), std::back_inserter(bead_pair_sep[bead_pair_id_2]), std::negate<double>());
                }
                else{
                    util->dist(bead->second, it->second, update_pair_sep[bead_pair_id_1], boxsize);
                    update_pair_sep[bead_pair_id_2].resize(0);
                    std::transform(update_pair_sep[bead_pair_id_1].begin(), update_pair_sep[bead_pair_id_1].end(), std::back_inserter(update_pair_sep[bead_pair_id_2]), std::negate<double>());
                }
            }
        }
    }
    
    const dVector& get_bead_separation(std::pair<size_t, size_t> bead_pair){
        if(update_pair_sep.find(bead_pair) != update_pair_sep.end())
            return update_pair_sep[bead_pair];
        else
            return bead_pair_sep[bead_pair];
    }
    
    void confirm_update(){
        
        for(boost::unordered_map<std::pair<size_t,size_t>, std::vector<double> >::iterator bead_sep = update_pair_sep.begin(); bead_sep!= update_pair_sep.end(); bead_sep++)
            bead_pair_sep[bead_sep->first] = bead_sep->second;
        
        update_pair_sep.clear();
    }
    
    void reject_update(){
        update_pair_sep.clear();
    }
    
    void set_update(bool status){
        update = status;
    }
    
    
private:
    
    std::vector<boost::unordered_map<size_t, std::vector<double>& > > beads_map;
    
    boost::unordered_map<std::pair<size_t,size_t>, std::vector<double> > bead_pair_sep;
    boost::unordered_map<std::pair<size_t,size_t>, std::vector<double> > update_pair_sep;
    
    bool update;
    
    boost::shared_ptr<Utility> util;
    
    double boxsize;
    
};

class Permutation_Separation_Table{
    
public:
    
    Permutation_Separation_Table(int multistep_dist, double boxsize = -1){
        util = boost::shared_ptr<Utility>(new Utility(0));
        this->multistep_dist = multistep_dist;
        this->boxsize = boxsize;
        update = true;
    };
    ~Permutation_Separation_Table(){};
    
    void add_bead(int slice, size_t key, std::vector<double>& location){
        if(slice >= beads_map.size())
            beads_map.resize(slice+1);
        beads_map[slice].insert(std::pair<size_t, std::vector<double>& >(key, location));
        std::pair<size_t, size_t> bead_pair = std::pair<size_t,size_t>(key,key);
        bead_pair_sep_vecs[bead_pair] = std::vector<double>(location.size(),0);
        bead_pair_sep[bead_pair] = 0.;
    }
    
    std::vector<double> get_bead(int slice, size_t key){
        return beads_map[slice][key];
    }
    
    void update_bead(int slice, size_t key){
        boost::unordered_map<size_t, std::vector<double>& >::iterator bead;
        bead = beads_map[slice].find(key);
        
        for(boost::unordered_map<size_t, std::vector<double>& >::iterator it = beads_map[(slice+multistep_dist)%beads_map.size()].begin(); it != beads_map[(slice+multistep_dist)%beads_map.size()].end(); it++){
            if(it->first != bead->first){
                std::pair<size_t, size_t> bead_pair_id_1 = std::pair<size_t,size_t>(bead->first,it->first);
                std::pair<size_t, size_t> bead_pair_id_2 = std::pair<size_t,size_t>(it->first,bead->first);
                
                if(update){
                    util->dist(bead->second, it->second, bead_pair_sep_vecs[bead_pair_id_1],boxsize);
                    bead_pair_sep_vecs[bead_pair_id_2] = bead_pair_sep_vecs[bead_pair_id_1];
                    bead_pair_sep[bead_pair_id_1] = inner_product(bead_pair_sep_vecs[bead_pair_id_1].begin(),bead_pair_sep_vecs[bead_pair_id_1].end(),bead_pair_sep_vecs[bead_pair_id_1].begin(),0.0);
                    bead_pair_sep[bead_pair_id_2] = bead_pair_sep[bead_pair_id_1];
                }
                else{
                    util->dist(bead->second, it->second, update_pair_sep_vecs[bead_pair_id_1],boxsize);
                    update_pair_sep_vecs[bead_pair_id_2] = update_pair_sep_vecs[bead_pair_id_1];
                    update_pair_sep[bead_pair_id_1] = inner_product(update_pair_sep_vecs[bead_pair_id_1].begin(),update_pair_sep_vecs[bead_pair_id_1].end(),update_pair_sep_vecs[bead_pair_id_1].begin(),0.0);
                    update_pair_sep[bead_pair_id_2] = update_pair_sep[bead_pair_id_1];

                }
            }
        }
        for(boost::unordered_map<size_t, std::vector<double>& >::iterator it = beads_map[(slice-multistep_dist+beads_map.size())%beads_map.size()].begin(); it != beads_map[(slice-multistep_dist+beads_map.size())%beads_map.size()].end(); it++){
            if(it->first != bead->first){
                std::pair<size_t, size_t> bead_pair_id_1 = std::pair<size_t,size_t>(bead->first,it->first);
                std::pair<size_t, size_t> bead_pair_id_2 = std::pair<size_t,size_t>(it->first,bead->first);
                
                if(update){
                    util->dist(bead->second, it->second, bead_pair_sep_vecs[bead_pair_id_1],boxsize);
                    bead_pair_sep_vecs[bead_pair_id_2] = bead_pair_sep_vecs[bead_pair_id_1];
                    bead_pair_sep[bead_pair_id_1] = inner_product(bead_pair_sep_vecs[bead_pair_id_1].begin(),bead_pair_sep_vecs[bead_pair_id_1].end(),bead_pair_sep_vecs[bead_pair_id_1].begin(),0.0);
                    bead_pair_sep[bead_pair_id_2] = bead_pair_sep[bead_pair_id_1];
                }
                else{
                    util->dist(bead->second, it->second, update_pair_sep_vecs[bead_pair_id_1],boxsize);
                    update_pair_sep_vecs[bead_pair_id_2] = update_pair_sep_vecs[bead_pair_id_1];
                    update_pair_sep[bead_pair_id_1] = inner_product(update_pair_sep_vecs[bead_pair_id_1].begin(),update_pair_sep_vecs[bead_pair_id_1].end(),update_pair_sep_vecs[bead_pair_id_1].begin(),0.0);
                    update_pair_sep[bead_pair_id_2] = update_pair_sep[bead_pair_id_1];
                    
                }
            }
        }
    }
    
    const dVector& get_perm_separation_vecs(std::pair<size_t, size_t> bead_pair){
        if(update_pair_sep_vecs.find(bead_pair) != update_pair_sep_vecs.end())
            return update_pair_sep_vecs[bead_pair];
        else
            return bead_pair_sep_vecs[bead_pair];
    }
    
    double get_perm_separation(std::pair<size_t, size_t> bead_pair){
        if(update_pair_sep.find(bead_pair) != update_pair_sep.end())
            return update_pair_sep[bead_pair];
        else
            return bead_pair_sep[bead_pair];
    }
    
    void confirm_update(){
        for(boost::unordered_map<std::pair<size_t,size_t>, std::vector<double> >::iterator bead_sep = update_pair_sep_vecs.begin(); bead_sep!= update_pair_sep_vecs.end(); bead_sep++)
            bead_pair_sep_vecs[bead_sep->first] = bead_sep->second;
        for(boost::unordered_map<std::pair<size_t,size_t>, double >::iterator bead_sep = update_pair_sep.begin(); bead_sep!= update_pair_sep.end(); bead_sep++)
            bead_pair_sep[bead_sep->first] = bead_sep->second;
        
        update_pair_sep_vecs.clear();
        update_pair_sep.clear();

    }
    
    void reject_update(){
        update_pair_sep_vecs.clear();
        update_pair_sep.clear();

    }
    
    void set_update(bool status){
        update = status;
    }
    
    
private:
    
    std::vector<boost::unordered_map<size_t, std::vector<double>& > > beads_map;
    
    boost::unordered_map<std::pair<size_t,size_t>, std::vector<double> > bead_pair_sep_vecs;
    boost::unordered_map<std::pair<size_t,size_t>, std::vector<double> > update_pair_sep_vecs;
    
    boost::unordered_map<std::pair<size_t,size_t>, double > bead_pair_sep;
    boost::unordered_map<std::pair<size_t,size_t>, double > update_pair_sep;

    
    bool update;
    int multistep_dist;
    
    boost::shared_ptr<Utility> util;
    
    double boxsize;
    
};

#endif /* lookuptable_hpp */
