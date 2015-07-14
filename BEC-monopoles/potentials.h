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
    potentials(){};
    ~potentials(){};
    double harmonicPotential(std::vector<double> loc, double m, double w);
    double free(){return 0;}
    double lj_int(double dist);
};



#endif /* defined(__PIMC__potentials__) */
