//
//  runner.hpp
//  PIMC-WORM
//
//  Created by Adith Ramamurti on 10/7/16.
//  Copyright © 2016 Adith Ramamurti. All rights reserved.
//

#ifndef runner_hpp
#define runner_hpp

#include "paths.hpp"
#include "parameters.hpp"
#include "moves.hpp"
#include "IO.hpp"
#include "paths.hpp"
#include "estimators.hpp"
#include <vector>
#include <boost/ptr_container/ptr_vector.hpp>
#include <stdio.h>

class Runner{
public:
    Runner(int &id, Parameters &params, MPI_Comm &local);
    ~Runner(){};
    void run(int &id, Parameters &params, IO &writer);
    void equilibrate(int &id, Parameters &params, Paths &paths, RNG &rng, Cos& cos);
    
private:
    boost::ptr_vector<Moves> moves;
    std::vector<std::vector<int> > move_choices;
    MPI_Comm local_comm;
    
};
#endif /* runner_hpp */
