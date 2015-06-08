//
//  paths.cpp
//  PIMCtest
//
//  Created by Adith Ramamurti on 5/20/15.
//  Copyright (c) 2015 Adith Ramamurti. All rights reserved.
//

#include "paths.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <sstream>


paths::paths(int procnum){
    ute = new utility(procnum);
    param = new parameters();
    //param->setT(0.5+procnum*0.5);
    last_chgd_part = {};
    multistep_dist = 16;
    printed = false;
    
    cout<< "Simulation Parameters:\nN      = \t" << param->getNumParticles() <<"\ntau    = \t" << param->gettau() << "\n" << "lambda =\t" << param->getlam() <<"\nT      = \t" << param->getT() << "\n\n";
    
    
    vector<vector<vector<double>>> beads(param->getNumTimeSlices(), vector<vector<double>>(param->getNumParticles(),vector<double>(param->getndim(), 0.0)));
    this->beads = beads;
    
    double offset[param->getNumParticles()][param->getndim()];
    for(int ptcl = 0; ptcl < param->getNumParticles(); ptcl++){
        nextConnection.push_back(ptcl);
        for(int ndim = 0; ndim < param->getndim(); ndim++){
            if(param->getBoxSize()!=-1)
                offset[ptcl][ndim] = param->getBoxSize()*(ute->randnormed(1)-0.5);
            else
                offset[ptcl][ndim] = ute->randnormed(1)-0.5;
        }
    }
    for(int slice = 0; slice < param->getNumTimeSlices(); slice++){
        for(int ptcl = 0; ptcl < param->getNumParticles(); ptcl++){
            for(int ndim = 0; ndim < param->getndim(); ndim++){
                this->beads[slice][ptcl][ndim]= offset[ptcl][ndim];//+0.5*(ute->randnormed(2)-1);
            }
        }
    }
    
    pot = new potentials();
    constPerms();
    
}

paths::~paths(){
    delete pot;
    delete ute;
    delete param;
}


void paths::constPerms(){
    vector<int> d(param->getNumParticles());
    int k;
    if(param->getNumParticles() <= 4)
        k = param->getNumParticles();
    else
        k = 4;
    
    vector<vector<int>> initPermList = {};
    iota(d.begin(),d.end(),0);
    
    do
    {
        vector<int> tempvec = {};
        for (int i = 0; i < k; i++)
        {
            tempvec.push_back(d[i]);
        }
        initPermList.push_back(tempvec);
        std::reverse(d.begin()+k,d.end());
    } while (next_permutation(d.begin(),d.end()));
    
    for(int i = 0; i < initPermList.size(); i++){
        vector<int> identity(param->getNumParticles());
        iota(identity.begin(),identity.end(),0);
        int displaced[k];
        for(int j = 0; j < k; j++){
            displaced[j] = initPermList[i][j];
        }
        sort(initPermList[i].begin(),initPermList[i].end());
        for(int j = 0; j < k; j++){
            identity[initPermList[i][j]] = displaced[j];
        }
        permList.push_back(identity);
    }
    
    sort(permList.begin(), permList.end());
    auto last = unique(permList.begin(), permList.end());
    permList.erase(last, permList.end());
    
    for(int i = 0; i < permList.size(); i++){
        vector<int> cp = {};
        for(int j = 0; j < permList[i].size(); j++)
            if(permList[i][j] != j)
                cp.push_back(j);
        sort(cp.begin(),cp.end());
        permPart.push_back(cp);
    }
    
    for(int slice = 0; slice < param->getNumTimeSlices(); slice ++){
        
        int end = (slice + multistep_dist)% param->getNumTimeSlices();
        vector<double> problist;
        
        vector<vector<double>> oldEndBeads(beads[end]);
        
        double oldAction = 0.0;
        for(int ptcl = 0; ptcl < param->getNumParticles(); ptcl++){
            oldAction += kineticAction(slice,multistep_dist);
        }
        
        for(int i = 0; i < permList.size(); i++){
            vector<int> oneperm = permList[i];
            int chdptcl = 0;
            for(int ptcl = 0; ptcl < param->getNumParticles(); ptcl++){
                beads[end][ptcl] = oldEndBeads[oneperm[ptcl]];
                if(ptcl != oneperm[ptcl]){
                    chdptcl += 1;
                }
            }
            double multvec[] = {1.0,20.0,100.0,400.0};
            
            double newAction = 0.0;
            for(int ptcl = 0; ptcl < param->getNumParticles(); ptcl++){
                newAction += kineticAction(slice,multistep_dist);
            }
            
            double multfac = 1.0;
            if(chdptcl != 0)
                multfac = multvec[chdptcl-1];
            problist.push_back(multfac*exp(-(newAction-oldAction)));
            beads[end] = oldEndBeads;
        }
        probList.push_back(problist);
    }
}


void paths::recompSingProb(vector<int> chdpart, int stslice){
    if(!chdpart.empty()){
    sort(chdpart.begin(),chdpart.end());
    vector<int> v(4);
    for(int perm = 0; perm < permList.size(); perm ++){
        vector<int>::iterator it = set_intersection(permPart[perm].begin(), permPart[perm].end(), chdpart.begin(),chdpart.end(),v.begin());
        if(it - v.begin() != 0){
            int end = (stslice + multistep_dist)% param->getNumTimeSlices();
            vector<vector<double>> oldEndBeads(beads[end]);
            
            double oldAction = 0.0;
            for(int ptcl = 0; ptcl < param->getNumParticles(); ptcl++){
                oldAction += kineticAction(stslice,multistep_dist);
            }
            double oldDens = exp(-oldAction);
            
            vector<int> oneperm = permList[perm];
            int chdptcl = 0;
            for(int ptcl = 0; ptcl < param->getNumParticles(); ptcl++){
                beads[end][ptcl] = oldEndBeads[oneperm[ptcl]];
                if(ptcl != oneperm[ptcl]){
                    chdptcl += 1;
                }
            }
            double multvec[] = {1.0,20.0,100.0,400.0};
            
            double action = 0.0;
            for(int ptcl = 0; ptcl < param->getNumParticles(); ptcl++){
                action += kineticAction(stslice,multistep_dist);
            }
            double newDens = exp(-action);
            double multfac = 1.0;
            if(chdptcl != 0)
                multfac = multvec[chdptcl-1];
            beads[end] = oldEndBeads;
            
            probList[stslice][perm]=multfac*(newDens/oldDens);
        }
    }
    }
}

double paths::vext(int slice, int ptcl){
    /*double vVal = 0.0;
    for(int i = 0; i < param->getNumParticles(); i++){
        if(i != ptcl){
            vector<double> distvec =ute->distance(beads[slice][ptcl], beads[slice][i], param->getBoxSize());

            vVal += pot->lj_int(sqrt(inner_product(distvec.begin(), distvec.end(),distvec.begin(), 0.0)));
        }
    }
    return vVal;*/
    return pot->harmonicPotential(beads[slice][ptcl], 1.0, 1.0);
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
        for(int ndim = 0; ndim < param->getndim(); ndim++){
            int slicep = (slice+dist);
            int ptcl2 = ptcl;
            if(slicep >= param->getNumTimeSlices()){
                ptcl2 = nextConnection[ptcl];
                slicep = slicep%param->getNumTimeSlices();
            }
            vector<double> distVec = ute->distance(beads[slicep][ptcl2], beads[slice][ptcl], -1);
            kin += 1/(2*param->getNumTimeSlices()*dist*param->gettau())*inner_product(distVec.begin(),distVec.end(),distVec.begin(),0.0);
        }
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
        int slicep1 = (slice+1);
        for(int ptcl = 0; ptcl < param->getNumParticles(); ptcl++){
            int ptcl2 = ptcl;
            if(slicep1 >= param->getNumTimeSlices()){
                ptcl2 = nextConnection[ptcl];
                slicep1 = slicep1%param->getNumTimeSlices();
            }
            vector<double> distVec = ute->distance(beads[slice][ptcl], beads[slicep1][ptcl2], -1);
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

void paths::print(){
    if(printed == false){
    ofstream outfile;
        for(int ptcl = 0; ptcl < param->getNumParticles(); ptcl++){
            stringstream sstm;
            sstm << "ptcl" << ptcl <<".csv";
            string result = sstm.str();
            outfile.open(result);
            for(int slice = 0; slice < param->getNumTimeSlices()*2; slice++){
                vector<double> a;
                if(slice < param->getNumTimeSlices())
                    a = beads[slice][ptcl];
                else
                    a = beads[slice%param->getNumTimeSlices()][nextConnection[ptcl]];
                for (int i = 0; i < a.size(); i++) {
                    stringstream strs;
                    outfile<<a[i];
                    string strel = strs.str();
                    strs.str("");
                    strs.clear();
                    outfile << strel+", ";
                }
                outfile<<"\n";
            }
            outfile.close();
        }
    }
    printed = true;
}
