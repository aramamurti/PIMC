//
//  utility.cpp
//  PIMC
//
//  Created by Adith Ramamurti on 5/20/15.
//  Copyright (c) 2015 Adith Ramamurti. All rights reserved.
//

#include "utility.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <numeric>
#include <cmath>


using namespace std;

string utility::currentDateTime() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    strftime(buf, sizeof(buf), "%Y-%m-%d--%H.%M.%S", &tstruct);
    return buf;
}

utility::utility(int procnum){
    gsl_rng_env_setup();
    
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    
    unsigned long int seed;
    seed = time(NULL)*procnum;
    gsl_rng_set(r,seed);

}

utility::~utility(){
    gsl_rng_free(r);
}

double utility::randnormed(int max){
    return gsl_rng_uniform (r) * max;
}

int utility::randint(int max){
    return (int)gsl_rng_uniform_int(r,max);
}

double utility::randgaussian(double width){
    return gsl_ran_gaussian_ziggurat(r, width);
}

double utility::vecavg(vector<double> v)
{
    double sum = accumulate(v.begin(), v.end(), 0.0);
    double mean = sum / v.size();
    
    return mean;
}

double utility::vecstd(vector<double> v){
    double mean = vecavg(v);
    vector<double> diff(v.size());
    std::transform(v.begin(), v.end(), diff.begin(),bind2nd(minus<double>(), mean));
    double sq_sum = inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
    double stdev = sqrt(sq_sum / v.size());
    
    return stdev;
}

vector<double> utility::location(vector<double> bead, double boxsize){
    int ndim = (int)bead.size();
    vector<double> loc;
    for(int i = 0; i < ndim; i++){
        if(boxsize != -1)
            loc.push_back(pbc(bead[i]+boxsize/2,boxsize)-boxsize/2);
        else
            loc.push_back(bead[i]);
    }
    return loc;
}
vector<double> utility::distance(vector<double> bead1, vector<double> bead2, double boxsize){
    int ndim = (int)bead1.size();
    vector<double> dist;
    bead1 = location(bead1,boxsize);
    bead2 = location(bead2,boxsize);
    for(int i = 0; i < ndim; i++){
        double dimdist = bead2[i]-bead1[i];
        if(boxsize != -1)
            dimdist = pbc(dimdist+boxsize/2, boxsize)-boxsize/2;
        dist.push_back(dimdist);
    }
    return dist;
}

double utility::pbc(double a, double b)
{
    if(b < 0)
        return fmod(-a, -b);
    double ret = fmod(a,b);
    if(ret < 0)
        ret+=b;
    return ret;
}
