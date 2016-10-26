//
//  utility.hpp
//
//
//  Created by Adith Ramamurti on 9/12/16.
//
//

#ifndef utility_h
#define utility_h

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <stdio.h>
#include <vector>
#include <deque>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iostream>



class RNG{
    
    const gsl_rng_type * T;
    gsl_rng * r;
    
public:
    
    RNG(){};
    
    void seed(int &id){
        gsl_rng_env_setup();
        T = gsl_rng_default;
        r = gsl_rng_alloc (T);
        unsigned long int seed;
        seed = time(NULL)*(id+1);
        gsl_rng_set(r,seed);
    }
    
    double randnormed(double max){
        return gsl_rng_uniform (r) * max;
    }
    
    int randint(int max){
        return (int)gsl_rng_uniform_int(r,max);
    }
    
    double randgaussian(double width){
        return gsl_ran_gaussian_ziggurat(r, width);
    }
};

extern double vector_avg(std::vector<double>& v);
extern double vector_std(std::vector<double>& v);
extern std::vector<double> vector_add(std::vector<double> a, std::vector<double> b);
extern double per_bound_cond(double a, double b);
extern void put_in_box(std::vector<double> &bead, double box_size);
extern void distance(const std::vector<double> &bead1, const std::vector<double> &bead2, std::vector<double> &dist, double box_size);
extern void average_loc(const std::vector<double>& bead1, const std::vector<double>& bead2, std::vector<double>& dstc, double box_size);
extern void average_loc_weighted(const std::vector<double>& bead1, const std::vector<double>& bead2, std::vector<double>& dstc, double box_size, int w1, int w2);


#endif /* utility_h */
