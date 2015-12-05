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

dVector Utility::location(dVector bead, double box_size){
    int ndim = (int)bead.size();
    dVector loc;
    for(int i = 0; i < ndim; i++){
        if(box_size != -1)
            loc.push_back(per_bound_cond(bead[i],box_size));
        else
            loc.push_back(bead[i]);
    }
    return loc;
}
dVector Utility::dist(ddVector beads, double box_size){
    int ndim = (int)beads[0].size();
    dVector dist;
    dVector bead1 = location(beads[0],box_size);
    dVector bead2 = location(beads[1],box_size);
    for(int i = 0; i < ndim; i++){
        double dimdist = bead2[i]-bead1[i];
        if(box_size != -1)
            dimdist = per_bound_cond(dimdist+box_size/2, box_size)-box_size/2;
        dist.push_back(dimdist);
    }
    return dist;
}

dVector Utility::avedist(ddVector beads, double box_size){
    dVector dstc = dist(beads, box_size);
    for(dVector::iterator it = dstc.begin(); it!= dstc.end(); it++){
        *it = *it/2 + beads[0][it-dstc.begin()];
    }
    return dstc;
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
