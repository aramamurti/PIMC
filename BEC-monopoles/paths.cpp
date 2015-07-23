//
//  paths.cpp
//  PIMCtest
//
//  Created by Adith Ramamurti on 5/20/15.
//  Copyright (c) 2015 Adith Ramamurti. All rights reserved.
//

#include "paths.h"


paths::paths(int procnum, std::ofstream &f)
:multvec{1.0,1.0,1.0,1.0}//{1.0,20.0,100.0,400.0}
{
    
    ute = new utility(procnum);
    param = new parameters();
    param->setT(1.0+procnum*0.15);
    multistep_dist = 8;
    numswap = 0;
    
    f<< "Simulation Parameters:\nN      = \t" << param->getNumParticles() <<"\nndim      = \t" << param->getndim() <<"\nBox Size      = \t" << param->getBoxSize() <<"\ntau    = \t" << param->gettau() << "\n" << "lambda =\t" << param->getlam() <<"\nT      = \t" << param->getT() << "\n" << std::endl;
    
    std::vector<std::vector<double>> offset(param->getNumParticles(), std::vector<double>(param->getndim(), 0.0));
    
    if(param->getBoxSize() == -1){
        for(int ptcl = 0; ptcl < param->getNumParticles(); ptcl++){
            for(int ndim = 0; ndim < param->getndim(); ndim++){
                offset[ptcl][ndim] = ute->randnormed(1)-0.5;
            }
        }
        beads = list_ptr(new LinkedList<std::vector<double>>());
        
        for(int slice = 0; slice < param->getNumTimeSlices(); slice++){
            for(int ptcl = 0; ptcl < param->getNumParticles(); ptcl++){
                beads->pushBack(offset[ptcl],ptcl);
            }
        }
        
    }
    else{
        offset.resize(0);
        double spacing = param->getBoxSize()/(pow(param->getNumParticles(),1/3.));
        int pps = (int)(pow(param->getNumParticles(),1/((double)param->getndim())));
        for(int i = 0; i < pps; i++){
            if(param->getndim() == 1){
                std::vector<double> pos;
                pos.push_back(i*spacing);
                offset.push_back(pos);
            }
            else{
                for(int j = 0; j < pps; j++){
                    if(param->getndim() == 1){
                        std::vector<double> pos;
                        pos.push_back(i*spacing);
                        pos.push_back(j*spacing);
                        offset.push_back(pos);
                    }
                    else{
                        for(int k = 0; k < pps; k++){
                            std::vector<double> pos;
                            pos.push_back(i*spacing);
                            pos.push_back(j*spacing);
                            pos.push_back(k*spacing);
                            offset.push_back(pos);
                        }
                    }
                }
            }
        }
        beads = list_ptr(new LinkedList<std::vector<double>>());
        
        for(int slice = 0; slice < param->getNumTimeSlices(); slice++){
            for(int ptcl = 0; ptcl < param->getNumParticles(); ptcl++){
                beads->pushBack(offset[ptcl],ptcl);
            }
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
            beads->setswap(identity,oneperm, slice, multistep_dist);
            beads->swap();
            double newAction = 0.0;
            for(int ptcl = 0; ptcl < param->getNumParticles(); ptcl++){
                newAction += kineticAction(slice,multistep_dist);
            }
            
            double multfac = 1.0;
            if(chdptcl != 0)
                multfac = multvec[chdptcl-1];
            permTot += exp(-(newAction - oldAction));
            problist.push_back(exp(-(newAction - oldAction)));
            
            beads->swap(true);
            
        }
        
        probList[slice]=problist;
        return permTot;
    }
    
    double paths::vext(int slice, int ptcl){
        
        double vVal = 0;
        if(param->getPots()[0])
            vVal += pot->harmonicPotential(beads->getOne(ptcl, slice), 1.0, 1.0);
        if(param->getPots()[1])
            for(int i = 0; i < param->getNumParticles(); i++){
                if(i != ptcl){
                    std::vector<double> distvec =ute->dist(beads->getPairSS(ptcl, i, slice), param->getBoxSize());
                    double dist = sqrt(inner_product(distvec.begin(), distvec.end(),distvec.begin(), 0.0));
                    vVal += pot->lj_int(dist);
                }
            }
        if(param->getPots()[2])
            for(int i = 0; i < param->getNumParticles(); i++){
                if(i != ptcl){
                    std::vector<double> distvec =ute->dist(beads->getPairSS(ptcl, i, slice), param->getBoxSize());
                    double dist = sqrt(inner_product(distvec.begin(), distvec.end(),distvec.begin(), 0.0));
                    vVal += pot->hardSphere(dist);
                }
            }
        
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
            std::vector<double> distVec = ute->dist(pair, param->getBoxSize());
            double ipdist =  inner_product(distVec.begin(),distVec.end(),distVec.begin(),0.0);
            kin += 1/(4*param->getlam()*param->gettau()*dist)*ipdist;
        }
        return kin;
    }
    
    /******************************************************
     Thermo Energy
     ******************************************************/
    
    double paths::potentialEnergy(){
        double PE = 0.0;
        for(int slice = 0; slice<param->getNumTimeSlices();slice++){
            for(int ptcl = 0; ptcl<param->getNumParticles(); ptcl++){
                if(param->getPots()[0])
                    PE += pot->harmonicPotential(beads->getOne(ptcl, slice), 1.0, 1.0);
                if(param->getPots()[1])
                    for(int i = ptcl+1; i < param->getNumParticles(); i++){
                        std::vector<double> distvec =ute->dist(beads->getPairSS(ptcl, i, slice), param->getBoxSize());
                        double dist = sqrt(inner_product(distvec.begin(), distvec.end(),distvec.begin(), 0.0));
                        PE += pot->lj_int(dist);
                    }
                if(param->getPots()[2])
                    for(int i = ptcl+1; i < param->getNumParticles(); i++){
                        std::vector<double> distvec =ute->dist(beads->getPairSS(ptcl, i, slice), param->getBoxSize());
                        double dist = sqrt(inner_product(distvec.begin(), distvec.end(),distvec.begin(), 0.0));
                        PE += pot->hardSphere(dist);
                    }
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
                std::vector<double> distVec = ute->dist(pair, param->getBoxSize());
                double dist = inner_product(distVec.begin(), distVec.end(), distVec.begin(), 0.0);
                tot -= norm*dist;
            }
        }
        
        double KE = 0.5*param->getndim()*param->getNumParticles()/param->gettau() +tot/param->getNumTimeSlices();
        
        return KE;
    }
    
    
    double paths::energy(){
        double energy = kineticEnergy()+potentialEnergy();
        return energy;
    }
    
    /*******************************************
     Virial Energy
     ********************************************/
    
    double paths::virialEnergy(){
        int virialWindow = param->getNumTimeSlices();
        
        double t1 = 3 *param->getndim() / (2*param->gettau()*virialWindow);
        
        double t2fac = 1 / (4 * param->getlam() * param->getNumTimeSlices() * virialWindow * pow(param->gettau(),2));
        double t2 = 0;
        for(int slice = 0; slice < param->getNumTimeSlices(); slice++){
            for(int ptcl = 0; ptcl < param->getNumParticles(); ptcl++){
                std::vector<double> iloc = beads->getOne(ptcl, slice);
                std::vector<std::vector<double>> p1;
                std::vector<std::vector<double>> p2;
                p1.push_back(iloc);
                p1.push_back(beads->getOne(ptcl, slice + virialWindow));
                
                p2.push_back(iloc);
                p2.push_back(beads->getOne(ptcl, slice + virialWindow-1));
                
                std::vector<double> gloc1 = ute->vecadd(iloc, ute->dist(p1, param->getBoxSize()));
                std::vector<double> gloc2 = ute->vecadd(iloc, ute->dist(p2, param->getBoxSize()));
                
                std::vector<std::vector<double>> p3;
                std::vector<std::vector<double>> p4;
                
                p3.push_back(iloc);
                p3.push_back(gloc1);
                
                p4.push_back(gloc1);
                p4.push_back(gloc2);
                
                
                std::vector<double> v1 = ute->dist(p3,param->getBoxSize());
                std::vector<double> v2 = ute->dist(p4,param->getBoxSize());
                t2 += inner_product(v1.begin(), v1.end(), v2.begin(), 0.0);
            }
        }
        t2*=t2fac;
        
        double t3fac = 1/(2*param->getNumTimeSlices());
        double t3 = 0;
        
        std::vector<std::vector<std::vector<double>>> rc;
        
        for(int slice = 0; slice < param->getNumTimeSlices(); slice++){
            std::vector<std::vector<double>> rcslice;
            for(int ptcl = 0; ptcl < param->getNumParticles(); ptcl++){
                std::vector<double> vsum(param->getndim(), 0.0);
                std::vector<double> iloc = beads->getOne(ptcl, slice);
                for(int gamma =0; gamma < virialWindow; gamma++){
                    std::vector<std::vector<double>> p1;
                    std::vector<std::vector<double>> p2;
                    p1.push_back(iloc);
                    p1.push_back(beads->getOne(ptcl, slice - gamma));
                    
                    p2.push_back(iloc);
                    p2.push_back(beads->getOne(ptcl, slice + gamma));
                    
                    std::vector<double> gloc1 = ute->vecadd(iloc, ute->dist(p1, param->getBoxSize()));
                    std::vector<double> gloc2 = ute->vecadd(iloc, ute->dist(p2, param->getBoxSize()));

                    ute->vecadd(vsum, ute->vecadd(gloc1, gloc2));
                }
                for(std::vector<double>::iterator it = vsum.begin(); it != vsum.end(); it++){
                    *it *= 1/(2.*virialWindow);
                }
                rcslice.push_back(vsum);
            }
            rc.push_back(rcslice);
        }
        
        if(param->getPots()[1])
            for(int slice = 0; slice < param->getNumTimeSlices(); slice++){
                for(int ptcl = 0; ptcl < param->getNumParticles(); ptcl++){
                    std::vector<double> gradSum(param->getndim(), 0.0);
                    for(int i = 0; i < param->getNumParticles(); i++){
                        if(i != ptcl){
                            std::vector<double> distvec =ute->dist(beads->getPairSS(i, ptcl, slice), param->getBoxSize());
                            double dist = sqrt(inner_product(distvec.begin(), distvec.end(),distvec.begin(), 0.0));
                            gradSum = ute->vecadd(gradSum, pot->grad_lj(distvec, dist));
                        }
                    }
                    std::vector<std::vector<double>> p1;
                    p1.push_back(beads->getOne(ptcl, slice));
                    p1.push_back(rc[slice][ptcl]);
                    std::vector<double> rcmr = ute->dist(p1,param->getBoxSize());
                    t3+= std::inner_product(rcmr.begin(), rcmr.end(), gradSum.begin(), 0.0);
                }
            }
        
        t3*= t3fac;
        
        double t4 = potentialEnergy();
        
        double energy = t1 + t2 + t3 + t4;
        return energy;
    }
    
    std::vector<int> paths::getCycles(){
        std::vector<int> cycles = beads->getCycles();
        return cycles;
    }
    
    void paths::putInBox(){
        for(int ptcl = 0; ptcl < param->getNumParticles(); ptcl++)
            for(int slice = 0; slice < param->getNumTimeSlices(); slice++){
                beads->setOne(ptcl, slice, ute->location(beads->getOne(ptcl, slice), param->getBoxSize()));
            }
        
    }
    
    
    
