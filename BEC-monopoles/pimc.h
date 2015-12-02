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
#include "estimators.hpp"

class PIMC{
public:
    PIMC(boost::shared_ptr<Path> path);
    ~PIMC(){};
    iVector run(int end_step, IO &writer, dVector &energytr, iiVector &cycleList);
    void set_up_moves(std::vector<bool> move_list);
    void set_up_estimators(std::vector<bool> estimator_list);
    void attempt_moves();
    void estimate();
    
private:
    boost::ptr_vector<Move_Base> moves;
    boost::ptr_vector<Estimator_Base> estimators;
    boost::shared_ptr<Path> path;
};


#endif
