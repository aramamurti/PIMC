//
//  utility.cpp
//
//
//  Created by Adith Ramamurti on 9/12/16.
//
//

#include "utility.hpp"

/*---------------------------------------------------------------------------------------------------*

This file contains the implementation of the various utility methods needed by the simulation.

*---------------------------------------------------------------------------------------------------*/

double vector_avg(std::vector<double>& v)
{
    double mean = 0;
    if(v.size()!=0){
        double sum = accumulate(v.begin(), v.end(), 0.0);
        mean = sum / v.size();
    }
    return mean;
}

double vector_std(std::vector<double>& v){
    double mean = vector_avg(v);
    double stdev = 0;
    if(v.size() != 0){
        std::vector<double> diff(v.size());
        std::transform(v.begin(), v.end(), diff.begin(),bind2nd(std::minus<double>(), mean));
        double sq_sum = inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
        stdev = sqrt(sq_sum / v.size());
    }
    return stdev;
}

std::vector<double> vector_add(std::vector<double> a, std::vector<double> b){
    std::vector<double> result;
    std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(result), [&](double l, double r)
                   {
                       return r+l;
                   });
    return result;
}

double per_bound_cond(double a, double b)
{
    if(b < 0)
        return fmod(-a, -b);
    double ret = fmod(a,b);
    if(ret < 0)
        ret+=b;
    return ret;
}

void put_in_box(std::vector<double> &bead, double box_size){
    if(box_size == -1)
        return;
    for(int i = 0; i < bead.size(); ++i){
        bead[i] = per_bound_cond(bead[i],box_size);
    }
}

void distance(const std::vector<double> &bead1, const std::vector<double> &bead2, std::vector<double> &dist, double box_size){
    dist.resize(bead1.size());
    for(int i = 0; i < dist.size(); ++i){
        double dimdist = bead2[i]-bead1[i];
        if(box_size != -1)
            dimdist = per_bound_cond(dimdist+box_size/2, box_size)-box_size/2;
        dist[i] = dimdist;
    }
}


void average_loc(const std::vector<double>& bead1, const std::vector<double>& bead2, std::vector<double>& dstc, double box_size){
    distance(bead1, bead2, dstc, box_size);
    for(std::vector<double>::iterator it = dstc.begin(); it!= dstc.end(); it++)
        *it = *it/2. + bead1[it-dstc.begin()];
}

void average_loc_weighted(const std::vector<double>& bead1, const std::vector<double>& bead2, std::vector<double>& dstc, double box_size, int w1, int w2){
    distance(bead1, bead2, dstc, box_size);
    for(std::vector<double>::iterator it = dstc.begin(); it!= dstc.end(); it++)
        *it = *it*(double(w1)/double(w1+w2)) + bead1[it-dstc.begin()];
}
