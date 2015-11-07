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


class moves{

public:
    moves(){};
    ~moves(){};
    
    bool comMove(Path* path, int ptcl);
    bool bisectionMoveHelper(Path* path, int ptcl);
    void bisectionMove(Path* path, int ptcl, int start, int m);
    vectori pickPermutation(Path* path, int start);
    
    
};

class Move_Base{
    
public:
    Move_Base(){};
    ~Move_Base(){};
};

class Center_of_Mass: public Move_Base{
    
public:
    Center_of_Mass(){};
    ~Center_of_Mass(){};
    
    
};

class Bisection: public Move_Base{
    
public:
    Bisection(){};
    ~Bisection(){};
    
};

class Insert: public Move_Base{
    
public:
    Insert(){};
    ~Insert(){};
    
};

class Remove: public Move_Base{
    
public:
    Remove(){};
    ~Remove(){};
    
};

class Open: public Move_Base{
    
public:
    Open(){};
    ~Open(){};
};

class Close: public Move_Base{
    
public:
    Close(){};
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
