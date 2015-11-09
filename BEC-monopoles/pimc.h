//
//  pimc.h
//  BEC-monopoles
//
//  Created by Adith Ramamurti on 6/8/15.
//  Copyright (c) 2015 Adith Ramamurti. All rights reserved.
//

#ifndef BEC_monopoles_pimc_h
#define BEC_monopoles_pimc_h

#include "moves.h"
#include "paths.h"
#include "uni_header.h"
#include <unistd.h>
#include "IO.hpp"

class Pimc{
public:
    Pimc();
    ~Pimc(){};
    std::vector<int> run(int end_step, Path* path, IO &writer, vectorf &energytr, vectorii &cycleList);
private:
    int numaccepts;
    int numacceptc;
    int numacceptb;
};


#endif
