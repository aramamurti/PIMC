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

    bool is_worm_move(){return worm_move;}
    bool is_worm_nec(){return worm_nec;}
    
    
protected:
    boost::shared_ptr<Path> path;
    
    int num_attempts;
    int num_accepts;
    
    int ptcl;
    int num_particles;
    
    double old_action;
    double new_action;
    
    bool worm_move;
    bool worm_nec;
    
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

class Insert: public Move_Base{
    
public:
    Insert(boost::shared_ptr<Path> path);
    ~Insert(){};
    
    Insert* clone() const{return new Insert(*this);}

    
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
    
    Remove* clone() const{return new Remove(*this);}

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
    
    Open* clone() const{return new Open(*this);}

    void attempt();
    void accept();
    bool check_move();
    
private:
    double mu_shift;
    int start_slice;
    int m;
    iVector rows_to_remove;
    
};

class Close: public Move_Base{
    
public:
    Close(boost::shared_ptr<Path> path);
    ~Close(){};
    
    Close* clone() const{return new Close(*this);}

    void attempt();
    void accept();
    void reject();
    bool check_move();
    
private:
    double mu_shift;
    int start_slice;
    int m;
    std::vector<std::pair<int, int> > ht;
    std::vector<std::pair<int, int> > start_end;


};

class Advance_Head: public Move_Base{
    
public:
    Advance_Head(boost::shared_ptr<Path> path);
    ~Advance_Head(){};
    
    Advance_Head* clone() const{return new Advance_Head(*this);}

    void attempt();
    void accept();
    void reject();
    void reject(int m);
    bool check_move();
    
private:
    double mu_shift;
    std::vector<std::pair<int, int> > start_end;

};

class Advance_Tail: public Move_Base{
    
public:
    Advance_Tail(boost::shared_ptr<Path> path);
    ~Advance_Tail(){};
    
    Advance_Tail* clone() const{return new Advance_Tail(*this);}

    void attempt();
    void accept();
    void reject();
    void reject(int m);
    bool check_move();
    
private:
    double mu_shift;
    std::vector<std::pair<int, int> > start_end;


};

class Recede_Head: public Move_Base{
    
public:
    Recede_Head(boost::shared_ptr<Path> path);
    ~Recede_Head(){};
    
    Recede_Head* clone() const{return new Recede_Head(*this);}

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
    
    Recede_Tail* clone() const{return new Recede_Tail(*this);}

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
    
    Swap_Head* clone() const{return new Swap_Head(*this);}

    void attempt();
    bool check_move();
    void reject();
    void accept();
    
private:
    int m;
    double sig_I;
    double sig_ksi;
    std::vector<std::pair<int, int> > start_end;
    iVector rows_to_remove;


};

class Swap_Tail: public Move_Base{
    
public:
    Swap_Tail(boost::shared_ptr<Path> path);
    ~Swap_Tail(){};
    
    Swap_Tail* clone() const{return new Swap_Tail(*this);}

    void attempt();
    bool check_move();
    void reject();
    void accept();
    
private:
    int m;
    double sig_I;
    double sig_ksi;
    std::vector<std::pair<int, int> > start_end;
    iVector rows_to_remove;
    
};

#endif /* defined(__PIMC__moves__) */
