//
//  constants.h
//  PIMC
//
//  Created by Adith Ramamurti on 5/14/15.
//  Copyright (c) 2015 Adith Ramamurti. All rights reserved.
//

#ifndef __PIMC__constants__
#define __PIMC__constants__


#endif /* defined(__PIMC__constants__) */
class constants{
    
private:
    double tau;
    double T;
    double hbar;
    double kb;
    double w;
    double m;
    double lambda;
    int ndim;
    
public:
    constants();
    ~constants();
    
    int getNdim(){return ndim;}
    double getHbar(){return hbar;}
    double getKb(){return kb;}
    double getOmega(){return w;}
    double getM(){return m;}    
};