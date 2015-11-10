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

class Move_Base{
    
public:
    Move_Base();
    Move_Base(boost::shared_ptr<Path> path);
    ~Move_Base(){};
    virtual Move_Base* clone() const = 0;
    
    void attempt();
    virtual void attempt(int ptcl) = 0;
    
    bool check_move();
    void accept();
    void reject();
    
    int get_num_accepts(){return num_accepts;}
    int get_num_attempts(){return num_attempts;}

    
    
protected:
    boost::shared_ptr<Path> path;
    int num_attempts;
    int num_accepts;
    
    float old_action;
    float new_action;

};

class Center_of_Mass: public Move_Base{
    
public:
    Center_of_Mass(boost::shared_ptr<Path> path);
    ~Center_of_Mass(){};
    
    Center_of_Mass* clone() const{return new Center_of_Mass(*this);}
    
    void attempt(int ptcl);
    void accept(int ptcl);
    
private:
    float delta;
    fVector shift;
    
};

class Bisection: public Move_Base{
    
public:
    Bisection(boost::shared_ptr<Path> path);
    ~Bisection(){};
    
    Bisection* clone() const{return new Bisection(*this);}
    
    void attempt(int ptcl);
    void accept(int ptcl);
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
    
    void attempt(int ptcl);
    void accept();
    void reject();
    iVector pick_perm(int start);
    
private:
    iVector permed_parts;
    
};

class Insert: public Move_Base{
    
public:
    Insert(boost::shared_ptr<Path> path) : Move_Base(path) {};
    ~Insert(){};
    
};

class Remove: public Move_Base{
    
public:
    Remove(boost::shared_ptr<Path> path) : Move_Base(path) {};
    ~Remove(){};
    
};

class Open: public Move_Base{
    
public:
    Open(boost::shared_ptr<Path> path) : Move_Base(path) {};
    ~Open(){};
};

class Close: public Move_Base{
    
public:
    Close(boost::shared_ptr<Path> path) : Move_Base(path) {};
    ~Close(){};
};

class Advance_Head: public Move_Base{
    
public:
    Advance_Head(){};
    ~Advance_Head(){};
};

class Advance_Tail: public Move_Base{
    
public:
    Advance_Tail(){};
    ~Advance_Tail(){};
};

class Recede_Head: public Move_Base{
    
public:
    Recede_Head(){};
    ~Recede_Head(){};
};

class Recede_Tail: public Move_Base{
    
public:
    Recede_Tail(){};
    ~Recede_Tail(){};
};

class Swap_Head: public Move_Base{
    
public:
    Swap_Head(){};
    ~Swap_Head(){};
};

class Swap_Tail: public Move_Base{
    
public:
    Swap_Tail(){};
    ~Swap_Tail(){};
};

#endif /* defined(__PIMC__moves__) */
