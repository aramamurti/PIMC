//
//  potentials.h
//  PIMC
//
//  Created by Adith Ramamurti on 5/14/15.
//  Copyright (c) 2015 Adith Ramamurti. All rights reserved.
//

#ifndef __PIMC__potentials__
#define __PIMC__potentials__

#include <stdio.h>
#include <vector>

using namespace std;

class potentials{
public:
    potentials(){};
    ~potentials(){};
    double harmonicPotential(double pos, double m, double w);
    double harmonicPotentialDeriv(vector<double> pos, double m, double w);
    double harmonicPotentialVir(vector<double> pos,double m, double w);
};



#endif /* defined(__PIMC__potentials__) */
