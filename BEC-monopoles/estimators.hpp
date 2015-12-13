//
//  estimators.hpp
//  BEC-monopoles
//
//  Created by Adith Ramamurti on 11/12/15.
//  Copyright © 2015 Adith Ramamurti. All rights reserved.
//

#ifndef estimators_hpp
#define estimators_hpp

#include <stdio.h>
#include "uni_header.h"
#include "path.h"
#include "potentials.h"
#include "permtable.hpp"

class Estimator_Base{
public:
    Estimator_Base(boost::shared_ptr<Path> path){
        this->path = path;
        utility = path->get_util();
        num_timeslices = path->get_parameters()->get_num_timeslices();
        ndim = path->get_parameters()->get_ndim();
    }
    ~Estimator_Base(){};
    virtual dVector estimate(){return dVector(0);}
    
protected:
    boost::shared_ptr<Path> path;
    boost::shared_ptr<Utility> utility;
    int num_timeslices;
    int ndim;

};

class Energy_Estimator: public Estimator_Base{
public:
    Energy_Estimator(boost::shared_ptr<Path> path, std::vector<bool> potentials): Estimator_Base(path){
        int i = 0;
        for(std::vector<bool>::iterator it = potentials.begin(); it != potentials.end(); it++){
            if(*it)
                switch(i){
                    case 0:
                        pot_funcs.push_back(boost::shared_ptr<Harmonic>(new Harmonic(1, 1)));
                        this->potentials.push_back(i);
                        break;
                    case 1:
                        pot_funcs.push_back(boost::shared_ptr<Lennard_Jones>(new Lennard_Jones()));
                        this->potentials.push_back(i);
                        break;
                    case 2:
                        pot_funcs.push_back(boost::shared_ptr<Aziz>(new Aziz()));
                        this->potentials.push_back(i);
                        break;
                    case 3:
                        pot_funcs.push_back(boost::shared_ptr<Coulomb>(new Coulomb(path->get_parameters()->get_num_particles(), pow(path->get_parameters()->get_box_size(),3),path->get_parameters()->get_coupling())));
                        this->potentials.push_back(i);
                        break;
                }
            i++;
        }
        
        norm = 1.0/(4.0*path->get_parameters()->get_lambda()*pow(path->get_parameters()->get_tau(),2));
        
    };
    ~Energy_Estimator(){};
    
    double potential_energy();
    double kinetic_energy();
    
    dVector estimate();
    
private:
    std::vector<boost::shared_ptr<Potential_Functions> > pot_funcs;
    std::vector<int> potentials;
    double norm;
};

class Permutation_Estimator: public Estimator_Base{
public:
    Permutation_Estimator(boost::shared_ptr<Path> path): Estimator_Base(path){};
    dVector estimate();
};

class Winding_Estimator: public Estimator_Base{
public:
    Winding_Estimator(boost::shared_ptr<Path> path): Estimator_Base(path){};
    dVector estimate();
};

#endif /* estimators_hpp */
