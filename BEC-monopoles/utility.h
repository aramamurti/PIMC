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


class Utility{
public:
    Utility(int procnum);
    ~Utility();
    std::string currentDateTime();
    double randnormed(int max);
    int randint(int max);
    double randgaussian(double width);
    double vecavg(dVector vec);
    double vecstd(dVector v);
    dVector vecadd(dVector a, dVector b);


    double per_bound_cond(double a, double b);
    int factorial(int n){return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;}
    int permutation(int n, int k){return factorial(n)/factorial(n-k);}
    dVector location(dVector bead, double box_size);
    dVector dist(ddVector beads, double box_size);
    dVector avedist(ddVector beads, double box_size);
    
private:
    const gsl_rng_type * T;
    gsl_rng * r;
};

#endif /* defined(__PIMC__utility__) */
