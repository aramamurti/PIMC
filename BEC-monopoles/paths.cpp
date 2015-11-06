//
//  paths.cpp
//  PIMCtest
//
//  Created by Adith Ramamurti on 5/20/15.
//  Copyright (c) 2015 Adith Ramamurti. All rights reserved.
//

#include "paths.h"

/****
 Constructor: sets up paths and parameters for each particle and constructs all objects necessary for the PIMC simulation steps.
 *****/
paths::paths(int procnum, std::ofstream &f)
:multvec{1.0,20.0,100.0, 400.0}{ // Multiplication factor for the swap probabilities
    
    //Declare all objects and variables, set all parameters that need to be modified from default (in parameters.h)
    ute = new utility(procnum);
    param = new parameters();
    param->setT(0.2+procnum*0.1);
    param->setNumTS((int)(40-20*param->getT()));
    multistep_dist = 8;
    last_start = 0;
    last_end = 0;
    pnum = procnum;
    pot = new potentials(param->getNumParticles(), pow(param->getBoxSize(),3));
    
    
    //Output to screen/file the parameters of the simulation
    f<< "Simulation Parameters:\nN      = \t" << param->getNumParticles() << "\nNumber of Time Slices   =\t" << param->getNumTimeSlices()<< "\nndim      = \t" << param->getndim() <<"\nBox Size      = \t" << param->getBoxSize() <<"\ntau    = \t" << param->gettau() << "\n" << "lambda =\t" << param->getlam() <<"\nT      = \t" << param->getT() << "\n" << std::endl;
    
    

    setup_beads();
    
    //If the particles are bosons, construct the permutation table
    if(param->isboson()){
        //Set up the probability table
        probList.resize(param->getNumTimeSlices());
        constr_perms(procnum);
    }
    
    }
    
    //Deconstructor
    paths::~paths(){
        delete pot;
        delete ute;
        delete param;
        beads.reset();
    }
    
    void paths::setup_beads(){
        //Set up the initial distribution of particles.
        vectorff offset(param->getNumParticles(), vectorf(param->getndim(), 0.0));
        
        if(param->getBoxSize() == -1){
            for(unsigned int ptcl = 0; ptcl < param->getNumParticles(); ptcl++){
                for(unsigned int ndim = 0; ndim < param->getndim(); ndim++){
                    offset[ptcl][ndim] = ute->randnormed(1)-0.5;
                }
            }
            
        }
        else{
            offset.resize(0);
            unsigned int pps = (int)ceil(pow(param->getNumParticles(),1/((float)param->getndim())));
            float spacing = param->getBoxSize()/pps;
            for(int i = 0; i < pps; i++){
                if(param->getndim() == 1){
                    vectorf pos;
                    pos.push_back(i*spacing);
                    offset.push_back(pos);
                }
                else{
                    for(int j = 0; j < pps; j++){
                        if(param->getndim() == 1){
                            vectorf pos;
                            pos.push_back(i*spacing);
                            pos.push_back(j*spacing);
                            offset.push_back(pos);
                        }
                        else{
                            for(int k = 0; k < pps; k++){
                                vectorf pos;
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
        
        beads = list_ptr(new LinkedList<vectorf>());
        
        for(int slice = 0; slice < param->getNumTimeSlices(); slice++){
            for(int ptcl = 0; ptcl < param->getNumParticles(); ptcl++){
                beads->pushBack(offset[ptcl],ptcl);
            }
        }
        
        if(param->ischarged()){
            int numcharges = param->getNumCharges();
            for(int i = 0; i < param->getNumParticles();i++){
                chgList.push_back(2*(i%numcharges)-1);
            }
        }
        
        //Make the beads a periodic chain
        beads->makeCircular();
        
    }
    
    //Construct permutation table
    void paths::constr_perms(int procnum){
        std::vector<int> d(param->getNumParticles());
        int k, maxPtcls = 3;
        if(param->getNumParticles() <= maxPtcls)
            k = param->getNumParticles();
        else
            k = maxPtcls;
        
        
        std::vector<std::vector<int>> initPermList(0);
        iota(d.begin(),d.end(),0);
        
        do
        {
            std::vector<int> tempvec(0);
            for (int i = 0; i < k; i++)
            {
                tempvec.push_back(d[i]);
            }
            initPermList.push_back(tempvec);
            std::reverse(d.begin()+k,d.end());
        } while (next_permutation(d.begin(),d.end()));
        
        for(std::vector<std::vector<int>>::iterator it = initPermList.begin(); it != initPermList.end(); it++){
            std::vector<int> identity(param->getNumParticles());
            iota(identity.begin(),identity.end(),0);
            int displaced[k];
            for(int j = 0; j < k; j++){
                displaced[j] = (*it)[j];
            }
            sort((*it).begin(),(*it).end());
            for(int j = 0; j < k; j++){
                identity[(*it)[j]] = displaced[j];
            }
            permList.push_back(identity);
        }
        
        sort(permList.begin(), permList.end());
        auto last = unique(permList.begin(), permList.end());
        permList.erase(last, permList.end());
        
        for(std::vector<std::vector<int>>::iterator it = permList.begin(); it != permList.end(); ){
            std::vector<int> cp(0);
            for(int j = 0; j < (*it).size(); j++)
                if((*it)[j] != j)
                    cp.push_back(j);
            sort(cp.begin(),cp.end());
            bool goodperm = true;
            if(param->ischarged() && param->getNumCharges()>1 && cp.size()>0){
                int chg1 = chgList[cp.front()];
                for(std::vector<int>::iterator it2 = cp.begin(); it2 != cp.end(); it2++){
                    if(chgList[*it2] != chg1){
                        goodperm = false;
                        break;
                    }
                }
            }
            if(goodperm){
                permPart.push_back(cp);
                ++it;
            }
            else{
                it = permList.erase(it);
            }
        }
        
        for(int ptcl = 0; ptcl < param->getNumParticles(); ptcl ++){
            std::vector<int> locs;
            
            for(int n = 0; n < permPart.size(); ++n)
            {
                auto i = find(permPart[n].begin(), permPart[n].end(), ptcl);
                if(permPart[n].end() != i)
                    locs.push_back(n);
            }
            permPartLoc.push_back(locs);
        }
        
        std::vector<int> ptcls(param->getNumParticles());
        iota(ptcls.begin(),ptcls.end(),0);
        
        
        for(int slice = 0; slice < param->getNumTimeSlices(); slice++){
            std::cout << procnum << ": Permutation table: slice " << slice <<std::endl;
            probList[slice].resize(permList.size());
            slice_perm_prob(ptcls, slice);
        }
    }
    
    
    float paths::slice_perm_prob(std::vector<int> ptcls, int slice){
        std::vector<int> permsToRecomp;
        for(std::vector<int>::iterator it = ptcls.begin(); it != ptcls.end(); it++){
            permsToRecomp.insert(permsToRecomp.end(), permPartLoc[*it].begin(), permPartLoc[*it].end());
        }
        std::vector<int>::iterator it;
        std::sort(permsToRecomp.begin(), permsToRecomp.end());
        it = std::unique (permsToRecomp.begin(), permsToRecomp.end());
        permsToRecomp.erase(it, permsToRecomp.end());
        
        
        float permTot = 0.0;
        std::vector<int> identity(param->getNumParticles());
        iota(identity.begin(),identity.end(),0);
        vectorf problist;
        
        float oldAction = 0.0;
        for(int ptcl = 0; ptcl < param->getNumParticles(); ptcl++){
            oldAction += kineticAction(slice,multistep_dist);
        }
        
        probList[slice][0] = 1;
        
        for(int j = 0; j < permsToRecomp.size(); j++){
            int i = permsToRecomp[j];
            std::vector<int> oneperm = permList[i];
            int chdptcl = 0;
            beads->setswap(identity,oneperm, slice, multistep_dist);
            beads->swap();
            
            chdptcl = (int)permPart[i].size();
            
            float newAction = 0.0;
            for(int ptcl = 0; ptcl < param->getNumParticles(); ptcl++){
                newAction += kineticAction(slice,multistep_dist);
            }
            
            float multfac = 1.0;
            if(chdptcl != 0)
                multfac = multvec[chdptcl-1];
            
            float kFac = multfac * exp(-(newAction - oldAction));
            probList[slice][i] = kFac;
            
            beads->swap(true);
        }
        
        for(vectorf::iterator it = probList[slice].begin(); it != probList[slice].end(); it++)
            permTot += *it;
        
        return permTot;
    }
    
    /********************
     External potential
     ******************/
    inline float paths::vext(int slice, int ptcl){
        
        float vVal = 0;
        if(param->getPots()[0])     //Harmonic
            vVal += pot->harmonicPotential(beads->getOne(ptcl, slice), 1.0, 1.0);
        if(param->getPots()[1]){    //LJ
            for(int i = 0; i < param->getNumParticles(); i++){
                if(i != ptcl){
                    vectorf distvec =ute->dist(beads->getPairSS(ptcl, i, slice), param->getBoxSize());
                    float dist = sqrt(inner_product(distvec.begin(), distvec.end(),distvec.begin(), 0.0));
                    vVal += pot->lj_int(dist);
                }
            }
            vVal = 0.5*vVal;
        }
        if(param->getPots()[2]){    //Hard Sphere
            for(int i = 0; i < param->getNumParticles(); i++){
                if(i != ptcl){
                    vectorf distvec =ute->dist(beads->getPairSS(ptcl, i, slice), param->getBoxSize());
                    float dist = sqrt(inner_product(distvec.begin(), distvec.end(),distvec.begin(), 0.0));
                    vVal += pot->hardSphere(dist);
                }
            }
            vVal = 0.5*vVal;
        }
        if(param->getPots()[3]){    //Aziz
            for(int i = 0; i < param->getNumParticles(); i++){
                if(i != ptcl){
                    vectorf distvec =ute->dist(beads->getPairSS(ptcl, i, slice), param->getBoxSize());
                    float dist = sqrt(inner_product(distvec.begin(), distvec.end(),distvec.begin(), 0.0));
                    vVal += pot->aziz_int(dist);
                }
            }
            vVal = 0.5*vVal;
        }
        if(param->getPots()[4]){    //Coulomb
            int chgi = chgList[ptcl];
            int nmax = pot->getnmax();
            for(int nx = -nmax; nx < nmax; nx++){
                for(int ny = -nmax; ny < nmax; ny++){
                    for(int nz = -nmax; nz < nmax; nz++){
                        for(int j = 0; j < param->getNumParticles(); j++){
                            
                        }
                    }
                }
            }
            
            vVal = 0.5*vVal;
        }
        
        return vVal;
    }
    
    float paths::potentialAction(int slice){
        float pot = 0;
        for(int ptcl = 0; ptcl < param->getNumParticles(); ptcl++){
            pot += vext(slice, ptcl);
        }
        return param->gettau()*pot;
    }
    
    float paths::kineticAction(int slice, int dist){
        float kin = 0;
        for(int ptcl = 0; ptcl < param->getNumParticles(); ptcl++){
            vectorff pair = beads->getPair(ptcl, slice, dist);
            vectorf distVec = ute->dist(pair, param->getBoxSize());
            float ipdist =  inner_product(distVec.begin(),distVec.end(),distVec.begin(),0.0);
            kin += 1/(4*param->getlam()*param->gettau()*dist)*ipdist;
        }
        return kin;
    }
    
    /******************************************************
     Thermo Energy
     ******************************************************/
    
    float paths::potentialEnergy(){
        float PE = 0.0;
        for(int slice = 0; slice<param->getNumTimeSlices();slice++){
            for(int ptcl = 0; ptcl<param->getNumParticles(); ptcl++){
                if(param->getPots()[0])
                    PE += pot->harmonicPotential(beads->getOne(ptcl, slice), 1.0, 1.0);
                if(param->getPots()[1])
                    for(int i = ptcl+1; i < param->getNumParticles(); i++){
                        vectorf distvec =ute->dist(beads->getPairSS(ptcl, i, slice), param->getBoxSize());
                        float dist = sqrt(inner_product(distvec.begin(), distvec.end(),distvec.begin(), 0.0));
                        PE += pot->lj_int(dist);
                    }
                if(param->getPots()[2])
                    for(int i = ptcl+1; i < param->getNumParticles(); i++){
                        vectorf distvec =ute->dist(beads->getPairSS(ptcl, i, slice), param->getBoxSize());
                        float dist = sqrt(inner_product(distvec.begin(), distvec.end(),distvec.begin(), 0.0));
                        PE += pot->hardSphere(dist);
                    }
                if(param->getPots()[3])
                    for(int i = ptcl+1; i < param->getNumParticles(); i++){
                        vectorf distvec =ute->dist(beads->getPairSS(ptcl, i, slice), param->getBoxSize());
                        float dist = sqrt(inner_product(distvec.begin(), distvec.end(),distvec.begin(), 0.0));
                        PE += pot->aziz_int(dist);
                    }
            }
        }
        PE = PE/param->getNumTimeSlices();
        return PE;
    }
    
    float paths::kineticEnergy(){
        float tot = 0.0;
        float norm = 1.0/(4.0*param->getlam()*pow(param->gettau(),2));
        for(int slice = 0; slice < param->getNumTimeSlices(); slice++){
            for(int ptcl = 0; ptcl < param->getNumParticles(); ptcl++){
                vectorff pair = beads->getPair(ptcl, slice, 1);
                vectorf distVec = ute->dist(pair, param->getBoxSize());
                float dist = inner_product(distVec.begin(), distVec.end(), distVec.begin(), 0.0);
                tot -= norm*dist;
            }
        }
        float KE = 0.5*param->getndim()*param->getNumParticles()/param->gettau() +tot/param->getNumTimeSlices();
        
        return KE;
    }
    
    float paths::energy(){
        float energy = kineticEnergy()+potentialEnergy();
        return energy;
    }
    
    std::vector<int> paths::getCycles(){
        std::vector<int> cycles = beads->getCycles();
        return cycles;
    }
    
    std::vector<int> paths::getWindingNumber(){
        vectorf dvectot(param->getndim(),0.0);
        for(int ptcl = 0; ptcl < param->getNumParticles(); ptcl++){
            for(int slice = 0; slice < param->getNumTimeSlices(); slice++){
                vectorff pair = beads->getPair(ptcl, slice, 1);
                vectorf distVec = ute->dist(pair, param->getBoxSize());
                dvectot = ute->vecadd(dvectot, distVec);
            }
        }
        std::vector<int> wnum(param->getndim(),0);
        for(vectorf::iterator it = dvectot.begin(); it != dvectot.end(); it++){
            wnum[it-dvectot.begin()] = (int)round(*it/param->getBoxSize());
        }
        return wnum;
    }
    
    void paths::putInBox(){
        if(param->getBoxSize() != -1)
            for(int ptcl = 0; ptcl < param->getNumParticles(); ptcl++)
                for(int slice = 0; slice < param->getNumTimeSlices(); slice++){
                    beads->setOne(ptcl, slice, ute->location(beads->getOne(ptcl, slice), param->getBoxSize()));
                }
        
    }
    
    
    
