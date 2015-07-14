//
//  pimc.h
//  BEC-monopoles
//
//  Created by Adith Ramamurti on 6/8/15.
//  Copyright (c) 2015 Adith Ramamurti. All rights reserved.
//

#ifndef BEC_monopoles_pimc_h
#define BEC_monopoles_pimc_h
#include "gnuplot_i.hpp"
class pimc{
public:
    pimc();
    ~pimc(){};
    std::vector<double> run(int numSteps, paths* path, Gnuplot &g);
private:
    int numaccepts;
    int numacceptc;
    int numacceptb;
};


#endif
