//
//  potentials.h
//  PIMC
//
//  Created by Adith Ramamurti on 5/14/15.
//  Copyright (c) 2015 Adith Ramamurti. All rights reserved.
//

#ifndef __PIMC__potentials__
#define __PIMC__potentials__

#include "uni_header.h"

class Potential_Functions{
public:
    Potential_Functions(){};
    Potential_Functions(int initna, double initvol){

    };
    
    ~Potential_Functions(){};
    
    
    virtual double potential_value(double){return 0;}
    virtual double potential_value(dVector){return 0;}
    virtual double potential_value(dVector, int, int, double){return 0;}
    
private:

};

class Harmonic: public Potential_Functions{
public:
    Harmonic(double m, double w){mass = m; omega = w;}
    ~Harmonic(){};
    
    double potential_value(dVector& loc){
        double pot_val = 0.0;
        for(int i = 0; i < loc.size(); i++){
            pot_val += 0.5*mass*pow(omega,2)*pow(loc[i],2);
        }
        return pot_val;
    }
    
private:
    double mass;
    double omega;
};

class Lennard_Jones: public Potential_Functions{
public:
    Lennard_Jones(){};
    ~Lennard_Jones(){};
    
     double potential_value(double dist){
        double eps = 10.22;
        double sig = 2.556;
        return 4*eps*(pow((sig/dist),12) - pow((sig/dist),6));
    }
};

class Aziz: public Potential_Functions{
public:
    Aziz(){};
    ~Aziz(){};
    
    inline double aziz_pcws(){
        if(dist >= 3.68335)
            return 1;
        else
            return exp(-pow((3.68335/dist-1),2));
    }
    
    double potential_value(double dist){
        this->dist = dist;
        double val = 10.8*(544850.4 * exp(-4.50018*dist)-(9424.94/pow(dist,10)+2556.63/pow(dist,8)+937.38/pow(dist,6))*aziz_pcws());
        return val;
    }
    
private:
    double dist;
    
};

class Coulomb: public Potential_Functions{
    
public:
    Coulomb(int initna, double initvol, double coupling){
        L = pow(initvol,1/3.);
        numatoms = initna;
        volume = initvol;
        alpha = sqrt(M_PI)*pow(numatoms/(pow(volume,2)),1/6.);
        coulcut = sqrt(p)/alpha;
        kcut = 2.*alpha*sqrt(p);
        nmax = floor(coulcut/L);
        kmax = ceil(kcut/(2.*M_PI*L));
        this->coupling = coupling;
    }
    
    ~Coulomb(){};
    
    void update_num_particles(int numatoms){
        this->numatoms = numatoms;
        alpha = sqrt(M_PI)*pow(numatoms/(pow(volume,2)),1/6.);
        coulcut = sqrt(p)/alpha;
        kcut = 2.*alpha*sqrt(p);
        nmax = floor(coulcut/L);
        kmax = ceil(kcut/(2.*M_PI*L));
    }
    
    inline double real_coulomb(double dist, int chgi, int chgj){
        double real = chgi*chgj * erfc(alpha*dist) / dist;
        return real;
    }
    
    inline double reci_coulomb(iVector kVec, dVector dVec, double sfac, double box_size, int chgi, int chgj){
        double kfac = 2*M_PI/box_size;
        double efac = 1;
        double k2 = 0;
        for(iVector::iterator it = kVec.begin(); it != kVec.end(); it++){
            efac = efac * cos(kfac * (*it) * dVec[it-kVec.begin()]);
            k2 += pow(kfac * *it,2);
        }
        double reci = sfac*(4*M_PI)/pow(box_size,3)*chgi*chgj*efac*exp(-k2*1/(4*alpha*alpha))/k2;
        return reci;
    }
    
    inline double self_coulomb(int chg, double box_size){
        return 2*alpha/sqrt(M_PI) * pow(chg,2);
    }
    
    double potential_value(dVector dist, int chgi, int chgj, double box_size){
        if(coupling == 0)
            return 0;
        
        double val = 0;
        for(int nx = -nmax; nx <= nmax; nx++)
            for(int ny = -nmax; ny <= nmax; ny++)
                for(int nz = -nmax; nz <= nmax; nz++){
                    double rSq = pow(dist[0] + nx*box_size, 2) + pow(dist[1] + ny*box_size, 2) + pow(dist[2] + nz*box_size, 2);
                    double r = sqrt(rSq);
                    if(rSq != 0){
                        val += real_coulomb(r, chgi, chgj);
                    }
                }
        
        for(int kx = 0; kx <= kmax; kx++)
            for(int ky = 0; ky <= kmax; ky++)
                for(int kz = 0; kz <= kmax; kz++){
                    iVector kVec(0);
                    kVec.push_back(kx);
                    kVec.push_back(ky);
                    kVec.push_back(kz);
                    double k2 = 0;
                    double kfac = 2*M_PI/box_size;
                    for(iVector::iterator it = kVec.begin(); it != kVec.end(); it++){
                        k2 += pow(kfac * *it,2);
                    }
                    if(k2 < pow(kcut,2) && (k2  != 0)){
                        int zeros = 0;
                        if(kx == 0)
                            zeros++;
                        if(ky == 0)
                            zeros++;
                        if(kz == 0)
                            zeros++;
                        
                        double sfac = 0.0;
                        switch(zeros){
                            case 1:
                                sfac = 8;
                                break;
                            case 2:
                                sfac = 4;
                                break;
                            case 3:
                                sfac = 2;
                                break;
                        }
                        val += reci_coulomb(kVec, dist, sfac, box_size, chgi, chgj);
                    }
                }
        
        if(dist[0] == 0 && dist[1] == 0 && dist[2] == 0)
            val -= self_coulomb(chgi, box_size);
        
        return coupling*val;
    }
    
    double getcoulcut(){return coulcut;}
    double getkcut(){return kcut;}
    int getnmax(){return nmax;}
    int getkmax(){return kmax;}


    
private:
    int numatoms;
    double volume;
    double ewaldError = 1E-16;
    double p = -log(ewaldError);
    double alpha;
    double coulcut;
    double kcut;
    int nmax;
    int kmax;
    double L;
    double coupling;
};

#endif /* defined(__PIMC__potentials__) */
