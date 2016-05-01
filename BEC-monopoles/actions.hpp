//
//  actions.hpp
//  BEC-monopoles
//
//  Created by Adith Ramamurti on 11/12/15.
//  Copyright © 2015 Adith Ramamurti. All rights reserved.
//

#ifndef actions_hpp
#define actions_hpp

#include <stdio.h>
#include "path.h"

class Action_Base{
public:
    Action_Base(boost::shared_ptr<Path> path){
        this->path = path;
        utility = path->get_util();
        num_timeslices = path->get_parameters()->get_num_timeslices();

    }
    ~Action_Base(){};
    double get_action(int slice, int dist);
protected:
    boost::shared_ptr<Path> path;
    boost::shared_ptr<Utility> utility;
    int num_timeslices;

};

class Potential_Action: public Action_Base{
public:
    Potential_Action(boost::shared_ptr<Path> path, std::vector<bool> potentials): Action_Base(path){
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
    }
    ~Potential_Action(){};
    
    double get_action(int slice, int dist, int start_omit = 0, int end_omit = 0);
    double potential_helper(int slice, int ptcl);
    double get_action_single_particle(int ptcl, int slice);
    double get_action_multiple_particles(iVector ptcls, int slice);
    
private:
    std::vector<boost::shared_ptr<Potential_Functions> > pot_funcs;
    std::vector<int> potentials;


};

#endif /* actions_hpp */
