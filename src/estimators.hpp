//
//  estimators.hpp
//  PIMC
//
//  Created by Adith Ramamurti on 10/25/16.
//  Copyright © 2016 Adith Ramamurti. All rights reserved.
//

#ifndef estimators_hpp
#define estimators_hpp

#include <stdio.h>
#include "paths.hpp"
#include "parameters.hpp"
#include "utility.hpp"

class Estimator{
public:
    Estimator(MPI_Comm& local){
        local_comm = local;
    };
    ~Estimator(){};
    void estimate(int& id, Paths& paths, Parameters &params, std::vector<double>& energy, std::vector<int>& winding, std::vector<int>& permutations);
private:
    MPI_Comm local_comm;
};


#endif /* estimators_hpp */
