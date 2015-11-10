//
//  utility.h
//  PIMC
//
//  Created by Adith Ramamurti on 5/20/15.
//  Copyright (c) 2015 Adith Ramamurti. All rights reserved.
//

#ifndef __PIMC__utility__
#define __PIMC__utility__

#include "uni_header.h"

#endif /* defined(__PIMC__utility__) */

class Utility{
public:
    Utility(int procnum);
    ~Utility();
    std::string currentDateTime();
    float randnormed(int max);
    int randint(int max);
    float randgaussian(float width);
    float vecavg(std::vector<float> vec);
    float vecstd(std::vector<float> v);
    std::vector<float> vecadd(std::vector<float> a, std::vector<float> b);


    float per_bound_cond(float a, float b);
    int factorial(int n){return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;}
    int permutation(int n, int k){return factorial(n)/factorial(n-k);}
    std::vector<float> location(std::vector<float> bead, float box_size);
    std::vector<float> dist(std::vector<std::vector<float> > beads, float box_size);
    std::vector<float> avedist(std::vector<std::vector<float> > beads, float box_size);
    
private:
    const gsl_rng_type * T;
    gsl_rng * r;
};