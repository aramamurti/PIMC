//
//  pimc.h
//  BEC-monopoles
//
//  Created by Adith Ramamurti on 6/8/15.
//  Copyright (c) 2015 Adith Ramamurti. All rights reserved.
//

#ifndef BEC_monopoles_pimc_h
#define BEC_monopoles_pimc_h

class pimc{
public:
    pimc();
    ~pimc(){};
    void run(int numSteps, paths* path, std::ofstream &f1, std::ofstream &f2, std::ofstream &f3, std::vector<double> &energytr, std::vector<std::vector<int>> &cycleList);
private:
    int numaccepts;
    int numacceptc;
    int numacceptb;
};


#endif
