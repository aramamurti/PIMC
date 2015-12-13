//
//  tables.cpp
//  BEC-monopoles
//
//  Created by Adith Ramamurti on 11/7/15.
//  Copyright © 2015 Adith Ramamurti. All rights reserved.
//

#include "permtable.hpp"

/****************************************************************************************
 
                            PERMUTATION TABLE CLASS
 
 ****************************************************************************************/

Permutation_Table::Permutation_Table(boost::shared_ptr<Path> path){
    this->path = path;
    ka = boost::shared_ptr<Kinetic_Action>(new Kinetic_Action(this->path));

    multistep_dist = path->get_multistep_dist();
    set_up_perms();
}

void Permutation_Table::set_up_perms(){
    prob_list.resize(path->get_parameters()->get_num_timeslices());
    
    iVector d(path->get_beads()->get_num_particles());
    int k, maxPtcls = 3;
    if(path->get_beads()->get_num_particles() <= maxPtcls)
        k = path->get_beads()->get_num_particles();
    else
        k = maxPtcls;
    
    iiVector init_perm_list(0);
    iota(d.begin(),d.end(),0);
    
    do
    {
        iVector tempvec(0);
        for (int i = 0; i < k; i++)
        {
            tempvec.push_back(d[i]);
        }
        init_perm_list.push_back(tempvec);
        std::reverse(d.begin()+k,d.end());
    } while (next_permutation(d.begin(),d.end()));
    
    for(iiVector::iterator it = init_perm_list.begin(); it != init_perm_list.end(); it++){
        iVector identity(path->get_parameters()->get_num_particles());
        iota(identity.begin(),identity.end(),0);
        int displaced[k];
        for(int j = 0; j < k; j++){
            displaced[j] = (*it)[j];
        }
        sort((*it).begin(),(*it).end());
        for(int j = 0; j < k; j++){
            identity[(*it)[j]] = displaced[j];
        }
        perm_list.push_back(identity);
    }
    
    sort(perm_list.begin(), perm_list.end());
    auto last = unique(perm_list.begin(), perm_list.end());
    perm_list.erase(last, perm_list.end());
    
    for(iiVector::iterator it = perm_list.begin(); it != perm_list.end(); ){
        iVector cp(0);
        for(int j = 0; j < (*it).size(); j++)
            if((*it)[j] != j)
                cp.push_back(j);
        sort(cp.begin(),cp.end());
        bool goodperm = true;
        if(path->get_parameters()->is_charged() && path->get_parameters()->get_num_chgs()>1 && cp.size()>0){
            int chg1 = (path->get_charge_list())[cp.front()];
            for(iVector::iterator it2 = cp.begin(); it2 != cp.end(); it2++){
                if((path->get_charge_list())[*it2] != chg1){
                    goodperm = false;
                    break;
                }
            }
        }
        if(goodperm){
            permed_parts.push_back(cp);
            ++it;
        }
        else{
            it = perm_list.erase(it);
        }
    }
    
    for(int ptcl = 0; ptcl < path->get_beads()->get_num_particles(); ptcl ++){
        iVector locs;
        
        for(int n = 0; n < permed_parts.size(); ++n)
        {
            auto i = find(permed_parts[n].begin(), permed_parts[n].end(), ptcl);
            if(permed_parts[n].end() != i)
                locs.push_back(n);
        }
        perm_part_loc.push_back(locs);
    }
    
    iVector ptcls(path->get_beads()->get_num_particles());
    iota(ptcls.begin(),ptcls.end(),0);
    
    
    for(int slice = 0; slice < path->get_parameters()->get_num_timeslices(); slice++){
        prob_list[slice].resize(perm_list.size());
        recalc_perms(ptcls, slice);
    }
}

double Permutation_Table::recalc_perms(iVector ptcls, int slice){
    iVector recomp_perms;
    
    for(iVector::iterator it = ptcls.begin(); it != ptcls.end(); it++){
        recomp_perms.insert(recomp_perms.end(), perm_part_loc[*it].begin(), perm_part_loc[*it].end());
    }
    
    iVector::iterator it;
    std::sort(recomp_perms.begin(), recomp_perms.end());
    it = std::unique (recomp_perms.begin(), recomp_perms.end());
    recomp_perms.erase(it, recomp_perms.end());
    
    
    double perm_tot = 0.0;
    iVector identity(path->get_beads()->get_num_particles());
    iota(identity.begin(),identity.end(),0);
    
    double old_action = 0.0;
    for(int ptcl = 0; ptcl < path->get_beads()->get_num_particles(); ptcl++){
        old_action += ka->get_action(slice,multistep_dist);
    }
    
    prob_list[slice][0] = 1;
    
    for(int j = 0; j < recomp_perms.size(); j++){
        int i = recomp_perms[j];
        iVector oneperm = perm_list[i];
        int chdptcl = 0;
        path->get_beads()->set_permutation(identity,oneperm, slice, multistep_dist);
        path->get_beads()->permute();
        
        chdptcl = (int)permed_parts[i].size();
        
        double new_action = 0.0;
        for(int ptcl = 0; ptcl < path->get_beads()->get_num_particles(); ptcl++){
            new_action += ka->get_action(slice,multistep_dist);
        }
        
        double multfac = 1.0;
        if(chdptcl != 0)
            multfac = path->multvec[chdptcl-1];
        
        double kFac = multfac * exp(-(new_action - old_action));
        prob_list[slice][i] = kFac;
        
        path->get_beads()->permute(true);
    }
    
    for(dVector::iterator it = prob_list[slice].begin(); it != prob_list[slice].end(); it++)
        perm_tot += *it;
    
    return perm_tot;
}

iVector Permutation_Table::pick_permutation(int ptcl, int start){
    dVector perm_weight = prob_list[start];
    iVector perm_loc = perm_part_loc[ptcl];
    
    double sum = perm_weight[0];
    for(iVector::iterator it = perm_loc.begin(); it != perm_loc.end(); it++)
        sum += perm_weight[*it];
    
    int choice = 0;
    double rn = path->get_util()->randnormed(sum);
    double probsum = 0.0;
    for(int i = 0; i < perm_loc.size()+1; i++){
        if(i == 0)
            probsum += perm_weight[0];
        else
            probsum += perm_weight[perm_loc[i-1]];
        if(rn<=probsum){
            if(i == 0)
                choice = 0;
            else
                choice = perm_loc[i-1];
            break;
        }
    }
    
    iVector chosenPerm = perm_list[choice];
    
    return chosenPerm;

}