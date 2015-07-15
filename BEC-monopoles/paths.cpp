//
//  paths.cpp
//  PIMCtest
//
//  Created by Adith Ramamurti on 5/20/15.
//  Copyright (c) 2015 Adith Ramamurti. All rights reserved.
//

#include "paths.h"


paths::paths(int procnum)
:multvec{1.0,1.0,1.0,1.0}//{1.0,20.0,100.0,400.0}
{
    printed = false;
    
    ute = new utility(procnum);
    param = new parameters();
    //param->setT(0.5+procnum*0.5);
    multistep_dist = 8;
    
    std::cout<< "Simulation Parameters:\nN      = \t" << param->getNumParticles() <<"\ntau    = \t" << param->gettau() << "\n" << "lambda =\t" << param->getlam() <<"\nT      = \t" << param->getT() << "\n\n";
    
    std::vector<std::vector<double>> offset(param->getNumParticles(), std::vector<double>(param->getndim(), 0.0));
    for(int ptcl = 0; ptcl < param->getNumParticles(); ptcl++){
        for(int ndim = 0; ndim < param->getndim(); ndim++){
            if(param->getBoxSize()!=-1){
                offset[ptcl][ndim] = param->getBoxSize()*(ute->randnormed(1)-0.5);
                std::cout<<offset[ptcl][ndim]<<std::endl;
            }
            else
                offset[ptcl][ndim] = ute->randnormed(1)-0.5;
        }
    }
    
    beads = list_ptr(new LinkedList<std::vector<double>>());
    
    for(int slice = 0; slice < param->getNumTimeSlices(); slice++){
        for(int ptcl = 0; ptcl < param->getNumParticles(); ptcl++){
            beads->pushBack(offset[ptcl],ptcl);
        }
    }
    
    probList.resize(param->getNumTimeSlices());
    
    beads->makeCircular();
        
    pot = new potentials();
    
    constPerms();
    }
    
    paths::~paths(){
        delete pot;
        delete ute;
        delete param;
        beads.reset();
    }
    
    
    void paths::constPerms(){
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
        
        for(std::vector<std::vector<int>>::iterator it = permList.begin(); it != permList.end(); it++){
            std::vector<int> cp(0);
            for(int j = 0; j < (*it).size(); j++)
                if((*it)[j] != j)
                    cp.push_back(j);
            sort(cp.begin(),cp.end());
            permPart.push_back(cp);
        }
    }
    
    
    double paths::recompSingProb(int slice){
        double permTot = 0.0;
        std::vector<int> identity(param->getNumParticles());
        iota(identity.begin(),identity.end(),0);
        std::vector<double> problist;
        
        double oldAction = 0.0;
        for(int ptcl = 0; ptcl < param->getNumParticles(); ptcl++){
            oldAction += kineticAction(slice,multistep_dist);
        }
        
        for(int i = 0; i < permList.size(); i++){
            std::vector<int> oneperm = permList[i];
            int chdptcl = 0;
            beads->swap(identity,oneperm, slice, multistep_dist);
            
            double newAction = 0.0;
            for(int ptcl = 0; ptcl < param->getNumParticles(); ptcl++){
                newAction += kineticAction(slice,multistep_dist);
            }

            double multfac = 1.0;
            if(chdptcl != 0)
                multfac = multvec[chdptcl-1];
            permTot += exp(-(newAction - oldAction));
            problist.push_back(exp(-(newAction - oldAction)));

            beads->reverseswap(identity,oneperm,slice,multistep_dist);
        }
        
        probList[slice]=problist;
        return permTot;
    }
    
    double paths::vext(int slice, int ptcl){
        
        double vVal = pot->harmonicPotential(beads->getOne(ptcl, slice), 1.0, 1.0);
//        if(param->getPots()[0])
//            vVal += pot->harmonicPotential(beads->getOne(ptcl, slice), 1.0, 1.0);
//        if(param->getPots()[1])
//            for(int i = 0; i < param->getNumParticles(); i++){
//                if(i != ptcl){
//                    std::vector<double> distvec =ute->distance(beads->getPair(ptcl, i, 0), param->getBoxSize());
//                    vVal += pot->lj_int(sqrt(inner_product(distvec.begin(), distvec.end(),distvec.begin(), 0.0)));
//                }
//            }
//        
//
        return vVal;

    }
    
    double paths::potentialAction(int slice){
        double pot = 0;
        for(int ptcl = 0; ptcl < param->getNumParticles(); ptcl++){
            pot += vext(slice, ptcl);
        }
        return param->gettau()*pot;
    }
    
    double paths::kineticAction(int slice, int dist){
        double kin = 0;
        for(int ptcl = 0; ptcl < param->getNumParticles(); ptcl++){
            std::vector<std::vector<double>> pair = beads->getPair(ptcl, slice, dist);
            std::vector<double> distVec = ute->dist(pair, -1);
            kin += 1/(4*param->getlam()*param->gettau()*dist)*inner_product(distVec.begin(),distVec.end(),distVec.begin(),0.0);
        }
        return kin;
    }
    
    double paths::potentialEnergy(){
        double PE = 0.0;
        for(int slice = 0; slice<param->getNumTimeSlices();slice++){
            for(int ptcl = 0; ptcl<param->getNumParticles(); ptcl++){
                PE += vext(slice, ptcl);
            }
        }
        PE = PE/param->getNumTimeSlices();
        return PE;
    }
    
    double paths::kineticEnergy(){
        double tot = 0.0;
        double norm = 1.0/(4.0*param->getlam()*pow(param->gettau(),2));
        for(int slice = 0; slice < param->getNumTimeSlices(); slice++){
            for(int ptcl = 0; ptcl < param->getNumParticles(); ptcl++){
                std::vector<std::vector<double>> pair = beads->getPair(ptcl, slice, 1);
                std::vector<double> distVec = ute->dist(pair, -1);
                double dist = inner_product(distVec.begin(), distVec.end(), distVec.begin(), 0.0);
                tot -= norm*dist;
            }
        }
        
        double KE = 0.5*param->getndim()*param->getNumParticles()/param->gettau() +tot/param->getNumTimeSlices();
        return KE;
    }
    
    
    double paths::energy(){
        double energy = kineticEnergy()+potentialEnergy();
        
        //std::cout<<kineticEnergy()<<std::endl;
        return energy;
    }
    
    void paths::setBeads(LinkedList<std::vector<double>> ll){
        beads->reset(ll);
    }
    
    
    
    
    
