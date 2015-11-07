//
//  pimc.h
//  BEC-monopoles
//
//  Created by Adith Ramamurti on 6/8/15.
//  Copyright (c) 2015 Adith Ramamurti. All rights reserved.
//

#ifndef BEC_monopoles_pimc_h
#define BEC_monopoles_pimc_h

class Pimc{
public:
    Pimc();
    ~Pimc(){};
    void run(int end_step, Path* path, std::ofstream &f1,std::ofstream &f2,std::ofstream &f3, vectorf &energytr, vectorii &cycleList);
private:
    int numaccepts;
    int numacceptc;
    int numacceptb;
};


#endif
