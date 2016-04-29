//
//  moves.h
//  PIMC
//
//  Created by Adith Ramamurti on 5/14/15.
//  Copyright (c) 2015 Adith Ramamurti. All rights reserved.
//

#ifndef __PIMC__moves__
#define __PIMC__moves__

#include "path.h"
#include "permtable.hpp"
#include "actions.hpp"

class Move_Base{
    
public:
    Move_Base(boost::shared_ptr<Path> path);
    virtual ~Move_Base(){};
    
    virtual Move_Base* clone() const = 0;
    
    virtual void attempt();
    
    virtual double get_delta(){return 0;}
    virtual void shift_delta(double shift){}
    
    void reset_acceptance_counters();
    
    bool check_move();
    void accept();
    void reject();
    
    int get_num_accepts(){return num_accepts;}
    int get_num_attempts(){return num_attempts;}
    
    std::string get_move_name(){return move_name;}

    bool is_perm_move(){return perm_move;}
    
    
protected:
    boost::shared_ptr<Path> path;
    
    int num_attempts;
    int num_accepts;
    
    int ptcl;
    int num_particles;
    
    double old_action;
    double new_action;
    
    bool perm_move;
    
    std::string move_name;
    
    iVector changed_particles;
    
    boost::shared_ptr<Potential_Action> pa;
    boost::shared_ptr<Kinetic_Action> ka;


};

class Center_of_Mass: public Move_Base{
    
public:
    Center_of_Mass(boost::shared_ptr<Path> path);
    ~Center_of_Mass(){};
    
    Center_of_Mass* clone() const{return new Center_of_Mass(*this);}
    
    void shift_delta(double shift);
    double get_delta(){return delta;}
    
    void attempt();
    void accept();
    
private:
    double delta;
    dVector shift;
    
};

class Bisection: public Move_Base{
    
public:
    Bisection(boost::shared_ptr<Path> path);
    ~Bisection(){};
    
    Bisection* clone() const{return new Bisection(*this);}
    
    void attempt();
    void accept();
    void reject();
    void level_move(int ptcl, int start, int m);
    
protected:
    int multistep_dist;
    int start;
    int end;
    
};

class Perm_Bisection: public Bisection{
public:
    Perm_Bisection(boost::shared_ptr<Path> path);
    ~Perm_Bisection(){};
    
    Perm_Bisection* clone() const{return new Perm_Bisection(*this);}
    
    void attempt();
    void accept();
    void reject();
    
private:
    iVector permed_parts;
    boost::shared_ptr<Permutation_Table> ptable;
    
};

#endif /* defined(__PIMC__moves__) */
