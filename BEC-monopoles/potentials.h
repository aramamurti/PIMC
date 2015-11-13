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
    Potential_Functions(int initna, float initvol){

    };
    
    ~Potential_Functions(){};
    
    
    virtual float potential_value(float){return 0;}
    virtual float potential_value(fVector){return 0;}
    
    virtual float real_coulomb(float, int, int){return 0;}
    virtual float reci_coulomb(fVector, int, float, int, int){return 0;}

private:

};

class Harmonic: public Potential_Functions{
public:
    Harmonic(float m, float w){mass = m; omega = w;}
    ~Harmonic(){};
    
    float potential_value(fVector loc){
        float pot_val = 0.0;
        for(int i = 0; i < loc.size(); i++){
            pot_val += 0.5*mass*pow(omega,2)*pow(loc[i],2);
        }
        return pot_val;
    }
    
private:
    float mass;
    float omega;
};

class Lennard_Jones: public Potential_Functions{
public:
    Lennard_Jones(){};
    ~Lennard_Jones(){};
    
     float potential_value(float dist){
        float eps = 10.22;
        float sig = 2.556;
        return 4*eps*(pow((sig/dist),12) - pow((sig/dist),6));
    }
};

class Aziz: public Potential_Functions{
public:
    Aziz(){};
    ~Aziz(){};
    
    inline float aziz_pcws(){
        if(dist >= 3.68335)
            return 1;
        else
            return exp(-pow((3.68335/dist-1),2));
    }
    
    float get_potential(float dist){
        this->dist = dist;
        float val = 10.8*(544850.4 * exp(-4.50018*dist)-(9424.94/pow(dist,10)+2556.63/pow(dist,8)+937.38/pow(dist,6))*aziz_pcws());
        return val;
    }
    
private:
    float dist;
    
};

class Coulomb: public Potential_Functions{
    
public:
    Coulomb(int initna, float initvol){
        L = pow(initvol,1/3.);
        numatoms = initna;
        volume = initvol;
        alpha = sqrt(M_PI)*pow(numatoms/(pow(volume,2)),1/6.);
        coulcut = sqrt(p)/alpha;
        kcut = 2.*alpha*sqrt(p);
        nmax = floor(coulcut/L);
        kmax = ceil(kcut/(2.*M_PI*L));
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
    
    inline float real_coulomb(float dist, int chgi, int chgj){
        float real = chgi*chgj * erfc(alpha*dist) / dist;
        return real;
    }
    
    inline float reci_coulomb(fVector kx, int sfac, float box_size, int chgi, int chgj){
        float kfac = 2* M_PI / box_size;
        float efac = 1;
        float k2 = 0;
        for(fVector::iterator it = kx.begin(); it != kx.end(); it++){
            efac = efac * cos(kfac * *it);
            k2 += pow(kfac * *it,2);
        }
        float reci = sfac/k2*exp(-k2*1/(4*alpha*alpha))*chgi *chgj*efac;
        return reci;
    }
    
    float getcoulcut(){return coulcut;}
    float getkcut(){return kcut;}
    int getnmax(){return nmax;}
    int getkmax(){return kmax;}


    
private:
    int numatoms;
    float volume;
    float ewaldError = 10E-7;
    float p = log(ewaldError);
    float alpha;
    float coulcut;
    float kcut;
    int nmax;
    int kmax;
    float L;
};

#endif /* defined(__PIMC__potentials__) */
