//
//  utility.h
//  PIMC
//
//  Created by Adith Ramamurti on 5/20/15.
//  Copyright (c) 2015 Adith Ramamurti. All rights reserved.
//

#ifndef __PIMC__utility__
#define __PIMC__utility__

#include <stdio.h>
#include <vector>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#endif /* defined(__PIMC__utility__) */

class utility{
public:
    utility(int procnum);
    ~utility();
    std::string currentDateTime();
    double randnormed(int max);
    int randint(int max);
    double randgaussian(double width);
    double vecavg(std::vector<double> vec);
    double vecstd(std::vector<double> v);
    double pbc(double a, double b);
    int factorial(int n){return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;}
    int permutation(int n, int k){return factorial(n)/factorial(n-k);}
    std::vector<double> location(std::vector<double> bead, double boxsize);
    std::vector<double> distance(std::vector<double> bead1, std::vector<double> bead2, double boxsize);
    
private:
    const gsl_rng_type * T;
    gsl_rng * r;
};