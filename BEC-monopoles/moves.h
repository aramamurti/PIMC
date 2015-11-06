//
//  moves.h
//  PIMC
//
//  Created by Adith Ramamurti on 5/14/15.
//  Copyright (c) 2015 Adith Ramamurti. All rights reserved.
//

#ifndef __PIMC__moves__
#define __PIMC__moves__

#include "paths.h"


class moves{

public:
    moves(){};
    ~moves(){};
    
    bool comMove(paths* path, int ptcl);
    bool bisectionMoveHelper(paths* path, int ptcl);
    void bisectionMove(paths* path, int ptcl, int start, int m);
    std::vector<int> pickPermutation(paths* path, int start);
    
    
};

#endif /* defined(__PIMC__moves__) */
