//
//  PIMC.h
//  BEC-monopoles
//
//  Created by Adith Ramamurti on 6/8/15.
//  Copyright (c) 2015 Adith Ramamurti. All rights reserved.
//

#ifndef BEC_monopoles_PIMC_h
#define BEC_monopoles_PIMC_h

#include "moves.h"
#include "paths.h"
#include "uni_header.h"
#include <unistd.h>
#include "IO.hpp"

class PIMC{
public:
    PIMC();
    ~PIMC(){};
    iVector run(int end_step, boost::shared_ptr<Path> path, IO &writer, fVector &energytr, iiVector &cycleList);
    void set_up_moves(boost::shared_ptr<Path> path, std::vector<bool> move_list);
    
private:
    boost::ptr_vector<Move_Base> moves;
};


#endif
