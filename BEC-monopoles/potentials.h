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

class potentials{
public:
    potentials(int initna, float initvol){
        L = pow(initvol,1/3.);
        numatoms = initna;
        volume = initvol;
        alpha = sqrt(M_PI)*pow(numatoms/(pow(volume,2)),1/6.);
        coulcut = sqrt(p)/alpha;
        kcut = 2.*alpha*sqrt(p);
        nmax = floor(coulcut/L);
        kmax = ceil(kcut/(2.*M_PI*L));
    };
    
    ~potentials(){};
    
    float getcoulcut(){return coulcut;}
    float getkcut(){return kcut;}
    int getnmax(){return nmax;}
    int getkmax(){return kmax;}
    
    float harmonicPotential(fVector loc, float m, float w);
    float free(){return 0;}
    float lj_int(float dist);
    float hardSphere(float dist);
    float aziz_int(float dist);
    float aziz_pcws(float dist);
    float real_coulomb(float dist, int chgi, int chgj);
    float reci_coulomb(fVector kx, int sfac, float box_size, int chgi, int chgj);

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
