//
//  pimc.h
//  PIMC
//
//  Created by Adith Ramamurti on 5/20/15.
//  Copyright (c) 2015 Adith Ramamurti. All rights reserved.
//

#ifndef __PIMC__pimc__
#define __PIMC__pimc__

#include <stdio.h>
#include <vector>
#include "paths.h"

using namespace std;

class pimc{
public:
    pimc();
    ~pimc(){};
    vector<double> run(int numSteps, paths* path, vector<bool> pmv);
    
    int numaccepts;
    int numacceptc;
};

#endif /* defined(__PIMC__pimc__) */
