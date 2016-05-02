//
//  utility.cpp
//  PIMC
//
//  Created by Adith Ramamurti on 5/20/15.
//  Copyright (c) 2015 Adith Ramamurti. All rights reserved.
//

#include "utility.h"


std::string Utility::currentDateTime() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    strftime(buf, sizeof(buf), "%Y-%m-%d--%H.%M.%S", &tstruct);
    return buf;
}

Utility::Utility(int procnum){
    gsl_rng_env_setup();
    
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    
    unsigned long int seed;
    seed = time(NULL)*(procnum+1);
    gsl_rng_set(r,seed);

}

Utility::~Utility(){
    gsl_rng_free(r);
}

double Utility::randnormed(double max){
    return gsl_rng_uniform (r) * max;
}

int Utility::randint(int max){
    return (int)gsl_rng_uniform_int(r,max);
}

double Utility::randgaussian(double width){
    return gsl_ran_gaussian_ziggurat(r, width);
}

double Utility::vecavg(dVector v)
{
    double mean = 0;
    if(v.size()!=0){
        double sum = accumulate(v.begin(), v.end(), 0.0);
        mean = sum / v.size();
    }
    return mean;
}

double Utility::vecstd(dVector v){
    double mean = vecavg(v);
    double stdev = 0;
    if(v.size() != 0){
        dVector diff(v.size());
        std::transform(v.begin(), v.end(), diff.begin(),bind2nd(std::minus<double>(), mean));
        double sq_sum = inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
        stdev = sqrt(sq_sum / v.size());
    }
    return stdev;
}

dVector Utility::vecadd(dVector a, dVector b){
    dVector result;
    std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(result), [&](double l, double r)
                   {
                       return r+l;
                   });
    
    return result;
}

void Utility::put_in_box(dVector &bead, double box_size){
    for(int i = 0; i < bead.size(); i++){
        if(box_size != -1)
            bead[i] = per_bound_cond(bead[i],box_size);
    }
}

void Utility::dist(dVector &bead1, dVector &bead2, dVector &dist, double box_size){
    int ndim = bead1.size();
    dist.resize(ndim);
    for(int i = 0; i < ndim; i++){
        double dimdist = bead2[i]-bead1[i];
        if(box_size != -1)
            dimdist = per_bound_cond(dimdist+box_size/2, box_size)-box_size/2;
        dist[i] = dimdist;
    }
}

void Utility::avedist(dVector& bead1, dVector& bead2, dVector& dstc, double box_size){
    dist(bead1, bead2, dstc, box_size);
    for(dVector::iterator it = dstc.begin(); it!= dstc.end(); it++){
        *it = *it/2 + bead1[it-dstc.begin()];
    }
}

double Utility::per_bound_cond(double a, double b)
{
    if(b < 0)
        return fmod(-a, -b);
    double ret = fmod(a,b);
    if(ret < 0)
        ret+=b;
    return ret;
}
