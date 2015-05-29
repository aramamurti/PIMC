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

using namespace std;

#endif /* defined(__PIMC__utility__) */

class utility{
public:
    utility(int procnum);
    ~utility();
    string currentDateTime();
    double randnormed(int max);
    int randint(int max);
    double randgaussian(double width);
    double vecavg(vector<double> vec);
    double vecstd(vector<double> v);
    void print(vector<double> vec);
    
private:
    const gsl_rng_type * T;
    gsl_rng * r;
};