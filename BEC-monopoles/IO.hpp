//
//  IO.hpp
//  BEC-monopoles
//
//  Created by Adith Ramamurti on 11/9/15.
//  Copyright © 2015 Adith Ramamurti. All rights reserved.
//

#ifndef IO_hpp
#define IO_hpp

#include "uni_header.h"
#include "parameters.h"


class IO{
public:
    IO(int world_rank){
        set_up_outfiles(world_rank);
    }
    
    void set_up_outfiles(int world_rank);
    void write_parameters(boost::shared_ptr<Parameters> params);
    void write_step_state(int step, fVector energy, iVector cycles, int num_particles, iVector wnum);
    void write_acceptance();
    void write_final(float energy, float energystd, int num_particles, iiVector cycles, iVector accept);
    void close(){    f1.close();
        f2.close();
        f3.close();
        f4.close();};
    
private:
    std::ofstream f1,f2,f3,f4;

};

#endif /* IO_hpp */
