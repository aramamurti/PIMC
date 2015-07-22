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
    std::vector<double> vecadd(std::vector<double> a, std::vector<double> b);


    double pbc(double a, double b);
    int factorial(int n){return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;}
    int permutation(int n, int k){return factorial(n)/factorial(n-k);}
    std::vector<double> location(std::vector<double> bead, double boxsize);
    std::vector<double> dist(std::vector<std::vector<double>> beads, double boxsize);
    std::vector<double> avedist(std::vector<std::vector<double>> beads, double boxsize);
    
private:
    const gsl_rng_type * T;
    gsl_rng * r;
};