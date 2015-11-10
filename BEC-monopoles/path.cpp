//
//  paths.cpp
//  PIMCtest
//
//  Created by Adith Ramamurti on 5/20/15.
//  Copyright (c) 2015 Adith Ramamurti. All rights reserved.
//

#include "path.h"

/****
 Constructor: sets up paths and parameters for each particle and constructs all objects necessary for the PIMC simulation steps.
 *****/
Path::Path(int procnum, IO &writer)
:multvec{1.0,20.0,100.0, 400.0}{ // Multiplication factor for the permute probabilities
    
    //Declare all objects and variables, set all parameters that need to be modified from default (in parameters.h)
    util = boost::shared_ptr<Utility>(new Utility(procnum));
    params = boost::shared_ptr<Parameters>(new Parameters());
//    params->set_T(0.2+procnum*0.1);
//    params->set_timeslices((int)(40-20*params->get_T()));
    
    writer.write_parameters(params);
    multistep_dist = 8;
    last_start = 0;
    last_end = 0;
    pnum = procnum;
    pot = boost::shared_ptr<potentials>(new potentials(params->get_num_particles(), pow(params->get_box_size(),3)));

    set_up_beads();
    
    //If the particles are bosons, construct the permutation table
    if(params->is_boson()){
        
        std::cout << procnum << ": Setting up paths and permutation table..." << std::endl;

        //Set up the probability table
        
        prob_list.resize(params->get_num_timeslices());
        constr_perms(procnum);
    }
    std::cout << procnum << ": Starting MC process ..." << std::endl;
    
    }
    
    //Deconstructor
    Path::~Path(){
        beads.reset();
    }
    
    void Path::set_up_beads(){
        //Set up the initial distribution of particles.
        ffVector offset(params->get_num_particles(), fVector(params->get_ndim(), 0.0));
        
        if(params->get_box_size() == -1){
            for(unsigned int ptcl = 0; ptcl < params->get_num_particles(); ptcl++){
                for(unsigned int ndim = 0; ndim < params->get_ndim(); ndim++){
                    offset[ptcl][ndim] = util->randnormed(1)-0.5;
                }
            }
            
        }
        else{
            offset.resize(0);
            unsigned int pps = (int)ceil(pow(params->get_num_particles(),1/((float)params->get_ndim())));
            float spacing = params->get_box_size()/pps;
            for(int i = 0; i < pps; i++){
                if(params->get_ndim() == 1){
                    fVector pos;
                    pos.push_back(i*spacing);
                    offset.push_back(pos);
                }
                else{
                    for(int j = 0; j < pps; j++){
                        if(params->get_ndim() == 1){
                            fVector pos;
                            pos.push_back(i*spacing);
                            pos.push_back(j*spacing);
                            offset.push_back(pos);
                        }
                        else{
                            for(int k = 0; k < pps; k++){
                                fVector pos;
                                pos.push_back(i*spacing);
                                pos.push_back(j*spacing);
                                pos.push_back(k*spacing);
                                offset.push_back(pos);
                            }
                        }
                    }
                }
            }
        }
        
        beads = list_ptr(new PathList<fVector>());
        
        for(int slice = 0; slice < params->get_num_timeslices(); slice++){
            for(int ptcl = 0; ptcl < params->get_num_particles(); ptcl++){
                beads->push_back(offset[ptcl],ptcl);
            }
        }
        
        if(params->is_charged()){
            int numcharges = params->get_num_chgs();
            for(int i = 0; i < params->get_num_particles();i++){
                charge_list.push_back(2*(i%numcharges)-1);
            }
        }
        
        //Make the beads a periodic chain
        beads->make_circular();
        
    }
    
    //Construct permutation table
    void Path::constr_perms(int procnum){
        iVector d(params->get_num_particles());
        int k, maxPtcls = 3;
        if(params->get_num_particles() <= maxPtcls)
            k = params->get_num_particles();
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
            iVector identity(params->get_num_particles());
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
            if(params->is_charged() && params->get_num_chgs()>1 && cp.size()>0){
                int chg1 = charge_list[cp.front()];
                for(iVector::iterator it2 = cp.begin(); it2 != cp.end(); it2++){
                    if(charge_list[*it2] != chg1){
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
        
        for(int ptcl = 0; ptcl < params->get_num_particles(); ptcl ++){
            iVector locs;
            
            for(int n = 0; n < permed_parts.size(); ++n)
            {
                auto i = find(permed_parts[n].begin(), permed_parts[n].end(), ptcl);
                if(permed_parts[n].end() != i)
                    locs.push_back(n);
            }
            perm_part_loc.push_back(locs);
        }
        
        iVector ptcls(params->get_num_particles());
        iota(ptcls.begin(),ptcls.end(),0);
        
        
        for(int slice = 0; slice < params->get_num_timeslices(); slice++){
            std::cout << procnum << ": Permutation table: slice " << slice <<std::endl;
            prob_list[slice].resize(perm_list.size());
            slice_perm_prob(ptcls, slice);
        }
    }
    
    
    float Path::slice_perm_prob(iVector ptcls, int slice){
        iVector recomp_perms;
        for(iVector::iterator it = ptcls.begin(); it != ptcls.end(); it++){
            recomp_perms.insert(recomp_perms.end(), perm_part_loc[*it].begin(), perm_part_loc[*it].end());
        }
        iVector::iterator it;
        std::sort(recomp_perms.begin(), recomp_perms.end());
        it = std::unique (recomp_perms.begin(), recomp_perms.end());
        recomp_perms.erase(it, recomp_perms.end());
        
        
        float permTot = 0.0;
        iVector identity(params->get_num_particles());
        iota(identity.begin(),identity.end(),0);

        float old_action = 0.0;
        for(int ptcl = 0; ptcl < params->get_num_particles(); ptcl++){
            old_action += kinetic_action(slice,multistep_dist);
        }
        
        prob_list[slice][0] = 1;
        
        for(int j = 0; j < recomp_perms.size(); j++){
            int i = recomp_perms[j];
            iVector oneperm = perm_list[i];
            int chdptcl = 0;
            beads->set_permutation(identity,oneperm, slice, multistep_dist);
            beads->permute();
            
            chdptcl = (int)permed_parts[i].size();
            
            float new_action = 0.0;
            for(int ptcl = 0; ptcl < params->get_num_particles(); ptcl++){
                new_action += kinetic_action(slice,multistep_dist);
            }
            
            float multfac = 1.0;
            if(chdptcl != 0)
                multfac = multvec[chdptcl-1];
            
            float kFac = multfac * exp(-(new_action - old_action));
            prob_list[slice][i] = kFac;
            
            beads->permute(true);
        }
        
        for(fVector::iterator it = prob_list[slice].begin(); it != prob_list[slice].end(); it++)
            permTot += *it;
        
        return permTot;
    }
    
    void Path::put_in_box(){
        if(params->get_box_size() != -1)
            for(int ptcl = 0; ptcl < params->get_num_particles(); ptcl++)
                for(int slice = 0; slice < params->get_num_timeslices(); slice++){
                    beads->set_bead_data(ptcl, slice, util->location(beads->get_bead_data(ptcl, slice), params->get_box_size()));
                }
        
    }

    
    /********************
     External potential
     ******************/
    inline float Path::vext(int slice, int ptcl){
        
        float vVal = 0;
        if(params->get_potentials()[0])     //Harmonic
            vVal += pot->harmonicPotential(beads->get_bead_data(ptcl, slice), 1.0, 1.0);
        if(params->get_potentials()[1]){    //LJ
            for(int i = 0; i < params->get_num_particles(); i++){
                if(i != ptcl){
                    fVector distvec =util->dist(beads->get_pair_same_slice(ptcl, i, slice), params->get_box_size());
                    float dist = sqrt(inner_product(distvec.begin(), distvec.end(),distvec.begin(), 0.0));
                    vVal += pot->lj_int(dist);
                }
            }
            vVal = 0.5*vVal;
        }
        if(params->get_potentials()[2]){    //Hard Sphere
            for(int i = 0; i < params->get_num_particles(); i++){
                if(i != ptcl){
                    fVector distvec =util->dist(beads->get_pair_same_slice(ptcl, i, slice), params->get_box_size());
                    float dist = sqrt(inner_product(distvec.begin(), distvec.end(),distvec.begin(), 0.0));
                    vVal += pot->hardSphere(dist);
                }
            }
            vVal = 0.5*vVal;
        }
        if(params->get_potentials()[3]){    //Aziz
            for(int i = 0; i < params->get_num_particles(); i++){
                if(i != ptcl){
                    fVector distvec =util->dist(beads->get_pair_same_slice(ptcl, i, slice), params->get_box_size());
                    float dist = sqrt(inner_product(distvec.begin(), distvec.end(),distvec.begin(), 0.0));
                    vVal += pot->aziz_int(dist);
                }
            }
            vVal = 0.5*vVal;
        }
        if(params->get_potentials()[4]){    //Coulomb
            int chgi = charge_list[ptcl];
            int nmax = pot->getnmax();
            for(int nx = -nmax; nx < nmax; nx++){
                for(int ny = -nmax; ny < nmax; ny++){
                    for(int nz = -nmax; nz < nmax; nz++){
                        for(int j = 0; j < params->get_num_particles(); j++){
                            
                        }
                    }
                }
            }
            
            vVal = 0.5*vVal;
        }
        
        return vVal;
    }
    
    float Path::potential_action(int slice){
        float pot = 0;
        for(int ptcl = 0; ptcl < params->get_num_particles(); ptcl++){
            pot += vext(slice, ptcl);
        }
        return params->get_tau()*pot;
    }
    
    float Path::kinetic_action(int slice, int dist){
        float kin = 0;
        for(int ptcl = 0; ptcl < params->get_num_particles(); ptcl++){
            ffVector pair = beads->get_pair_same_path(ptcl, slice, dist);
            fVector distVec = util->dist(pair, params->get_box_size());
            float ipdist =  inner_product(distVec.begin(),distVec.end(),distVec.begin(),0.0);
            kin += 1/(4*params->get_lambda()*params->get_tau()*dist)*ipdist;
        }
        return kin;
    }
    
    /************************************************************************************
     Estimators: Energy, Permutation, Winding Number
     *************************************************************************************/
    
    float Path::potential_energy(){
        float PE = 0.0;
        for(int slice = 0; slice<params->get_num_timeslices();slice++){
            for(int ptcl = 0; ptcl<params->get_num_particles(); ptcl++){
                if(params->get_potentials()[0])
                    PE += pot->harmonicPotential(beads->get_bead_data(ptcl, slice), 1.0, 1.0);
                if(params->get_potentials()[1])
                    for(int i = ptcl+1; i < params->get_num_particles(); i++){
                        fVector distvec =util->dist(beads->get_pair_same_slice(ptcl, i, slice), params->get_box_size());
                        float dist = sqrt(inner_product(distvec.begin(), distvec.end(),distvec.begin(), 0.0));
                        PE += pot->lj_int(dist);
                    }
                if(params->get_potentials()[2])
                    for(int i = ptcl+1; i < params->get_num_particles(); i++){
                        fVector distvec =util->dist(beads->get_pair_same_slice(ptcl, i, slice), params->get_box_size());
                        float dist = sqrt(inner_product(distvec.begin(), distvec.end(),distvec.begin(), 0.0));
                        PE += pot->hardSphere(dist);
                    }
                if(params->get_potentials()[3])
                    for(int i = ptcl+1; i < params->get_num_particles(); i++){
                        fVector distvec =util->dist(beads->get_pair_same_slice(ptcl, i, slice), params->get_box_size());
                        float dist = sqrt(inner_product(distvec.begin(), distvec.end(),distvec.begin(), 0.0));
                        PE += pot->aziz_int(dist);
                    }
            }
        }
        PE = PE/params->get_num_timeslices();
        return PE;
    }
    
    float Path::kinetic_energy(){
        float tot = 0.0;
        float norm = 1.0/(4.0*params->get_lambda()*pow(params->get_tau(),2));
        for(int slice = 0; slice < params->get_num_timeslices(); slice++){
            for(int ptcl = 0; ptcl < params->get_num_particles(); ptcl++){
                ffVector pair = beads->get_pair_same_path(ptcl, slice, 1);
                fVector distVec = util->dist(pair, params->get_box_size());
                float dist = inner_product(distVec.begin(), distVec.end(), distVec.begin(), 0.0);
                tot -= norm*dist;
            }
        }
        float KE = 0.5*params->get_ndim()*params->get_num_particles()/params->get_tau() +tot/params->get_num_timeslices();
        
        return KE;
    }
    
    float Path::energy(){
        float energy = kinetic_energy()+potential_energy();
        return energy;
    }
    
    iVector Path::get_cycles(){
        iVector cycles = beads->get_cycles();
        return cycles;
    }
    
    iVector Path::get_winding_number(){
        fVector dvectot(params->get_ndim(),0.0);
        for(int ptcl = 0; ptcl < params->get_num_particles(); ptcl++){
            for(int slice = 0; slice < params->get_num_timeslices(); slice++){
                ffVector pair = beads->get_pair_same_path(ptcl, slice, 1);
                fVector distVec = util->dist(pair, params->get_box_size());
                dvectot = util->vecadd(dvectot, distVec);
            }
        }
        iVector wnum(params->get_ndim(),0);
        for(fVector::iterator it = dvectot.begin(); it != dvectot.end(); it++){
            wnum[it-dvectot.begin()] = (int)round(*it/params->get_box_size());
        }
        return wnum;
    }
    
    
