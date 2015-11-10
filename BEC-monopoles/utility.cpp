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

float Utility::randnormed(int max){
    return gsl_rng_uniform (r) * max;
}

int Utility::randint(int max){
    return (int)gsl_rng_uniform_int(r,max);
}

float Utility::randgaussian(float width){
    return gsl_ran_gaussian_ziggurat(r, width);
}

float Utility::vecavg(fVector v)
{
    float sum = accumulate(v.begin(), v.end(), 0.0);
    float mean = sum / v.size();
    
    return mean;
}

float Utility::vecstd(fVector v){
    float mean = vecavg(v);
    fVector diff(v.size());
    std::transform(v.begin(), v.end(), diff.begin(),bind2nd(std::minus<float>(), mean));
    float sq_sum = inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
    float stdev = sqrt(sq_sum / v.size());
    
    return stdev;
}

fVector Utility::vecadd(fVector a, fVector b){
    fVector result;
    std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(result), [&](float l, float r)
                   {
                       return r+l;
                   });
    
    return result;
}

fVector Utility::location(fVector bead, float box_size){
    int ndim = (int)bead.size();
    fVector loc;
    for(int i = 0; i < ndim; i++){
        if(box_size != -1)
            loc.push_back(per_bound_cond(bead[i],box_size));
        else
            loc.push_back(bead[i]);
    }
    return loc;
}
fVector Utility::dist(ffVector beads, float box_size){
    int ndim = (int)beads[0].size();
    fVector dist;
    fVector bead1 = location(beads[0],box_size);
    fVector bead2 = location(beads[1],box_size);
    for(int i = 0; i < ndim; i++){
        float dimdist = bead2[i]-bead1[i];
        if(box_size != -1)
            dimdist = per_bound_cond(dimdist+box_size/2, box_size)-box_size/2;
        dist.push_back(dimdist);
    }
    return dist;
}

fVector Utility::avedist(ffVector beads, float box_size){
    fVector dstc = dist(beads, box_size);
    for(fVector::iterator it = dstc.begin(); it!= dstc.end(); it++){
        *it = *it/2 + beads[0][it-dstc.begin()];
    }
    return dstc;
}

float Utility::per_bound_cond(float a, float b)
{
    if(b < 0)
        return fmod(-a, -b);
    float ret = fmod(a,b);
    if(ret < 0)
        ret+=b;
    return ret;
}
