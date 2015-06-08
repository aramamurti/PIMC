//
//  moves.h
//  PIMC
//
//  Created by Adith Ramamurti on 5/14/15.
//  Copyright (c) 2015 Adith Ramamurti. All rights reserved.
//

#ifndef __PIMC__moves__
#define __PIMC__moves__

#include <stdio.h>
#include <vector>
#include "paths.h"

class moves{

public:
    moves(){};
    ~moves(){};
    
    bool comMove(paths* path, int ptcl);
    bool stagingMoveHelper(paths* path, int ptcl);
    void stagingMove(paths* path, int ptcl, int alpha_start, int alpha_end, int m);
    void bisectionMove(paths* path, int ptcl, int alpha_start, int alpha_end, int m);
    bool bisectionMoveHelper(paths* path, int ptcl);
    std::vector<int> pickPermutation(paths* path, int alpha_start, int alpha_end);
};

#endif /* defined(__PIMC__moves__) */
