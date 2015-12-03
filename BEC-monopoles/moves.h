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
    ~Move_Base(){};
    
    virtual Move_Base* clone() const = 0;
    
    virtual void attempt();
    
    bool check_move();
    void accept();
    void reject();
    
    int get_num_accepts(){return num_accepts;}
    int get_num_attempts(){return num_attempts;}
    
    std::string get_move_name(){return move_name;}

    
    
protected:
    boost::shared_ptr<Path> path;
    
    int num_attempts;
    int num_accepts;
    
    int ptcl;
    
    double old_action;
    double new_action;
    
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

class Insert: public Move_Base{
    
public:
    Insert(boost::shared_ptr<Path> path);
    ~Insert(){};
    void attempt();
    void accept();
    void reject();
    bool check_move();
    
private:
    double mu_shift;

    
};

class Remove: public Move_Base{
    
public:
    Remove(boost::shared_ptr<Path> path);
    ~Remove(){};
    void attempt();
    void accept();
    void reject();
    bool check_move();

private:
    double mu_shift;
    
};

class Open: public Move_Base{
    
public:
    Open(boost::shared_ptr<Path> path);
    ~Open(){};
    void attempt();
    void accept();
    bool check_move();
    
private:
    double mu_shift;
    int start_slice;
    int m;
    
};

class Close: public Move_Base{
    
public:
    Close(boost::shared_ptr<Path> path);
    ~Close(){};
    void attempt();
    void accept();
    void reject();
    bool check_move();
    
private:
    double mu_shift;
    int start_slice;
    int m;
    std::vector<std::pair<int, int> > ht;
};

class Advance_Head: public Move_Base{
    
public:
    Advance_Head(boost::shared_ptr<Path> path);
    ~Advance_Head(){};
    void attempt();
    void accept();
    void reject();
    void reject(int m);
    bool check_move();
    
private:
    double mu_shift;

};

class Advance_Tail: public Move_Base{
    
public:
    Advance_Tail(boost::shared_ptr<Path> path);
    ~Advance_Tail(){};
    void attempt();
    void accept();
    void reject();
    void reject(int m);
    bool check_move();
    
private:
    double mu_shift;


};

class Recede_Head: public Move_Base{
    
public:
    Recede_Head(boost::shared_ptr<Path> path);
    ~Recede_Head(){};
    void attempt();
    void accept(int m);
    void reject();
    bool check_move();
    
private:
    double mu_shift;

};

class Recede_Tail: public Move_Base{
    
public:
    Recede_Tail(boost::shared_ptr<Path> path);
    ~Recede_Tail(){};
    void attempt();
    void accept(int m);
    void reject();
    bool check_move();
    
private:
    double mu_shift;

};

class Swap_Head: public Move_Base{
    
public:
    Swap_Head(boost::shared_ptr<Path> path);
    ~Swap_Head(){};
};

class Swap_Tail: public Move_Base{
    
public:
    Swap_Tail(boost::shared_ptr<Path> path);
    ~Swap_Tail(){};
};

#endif /* defined(__PIMC__moves__) */
