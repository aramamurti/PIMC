//
//  IO.hpp
//  BEC-monopoles
//
//  Created by Adith Ramamurti on 11/9/15.
//  Copyright © 2015 Adith Ramamurti. All rights reserved.
//

#ifndef IO_hpp
#define IO_hpp

#include "parameters.hpp"
#include "moves.hpp"
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <vector>




class IO{
public:
    IO(int& id, int world_rank, MPI_Comm &local){
        local_comm = local;
        if(id == 0)
            set_up_outfiles(world_rank);
    }
    ~IO(){}
    
    void set_up_outfiles(int& id);
    void write_parameters(Parameters &params);
    void write_equil_parameters(Parameters& params, double delta);
    void write_step_state(int step, int num_particles, std::vector<double>& energy, std::vector<int>& cycles,std::vector<int>& wnum);
    void write_final(int counter, double particles, double energy_th, double energy_th_std, double energy_v, double energy_v_std, std::vector<std::vector<int> >& cycles, boost::ptr_vector<Moves> &moves);
    void write_acceptances(int counter,boost::ptr_vector<Moves> &moves);
    void path_dump(int &id, int count, Paths& paths, Parameters &params);
    void close(){
        f1.close();
        f2.close();
        f3.close();
        f4.close();
        f5.close();
        pathf << "}" << std::endl;
        pathf.close();
        pathf2 << "}" << std::endl;
        pathf2.close();
    };
    
    void read_parameters(std::string infile, Parameters& params);
    
private:
    std::ofstream f1,f2,f3,f4,f5;
    std::ofstream pathf, pathf2;
    MPI_Comm local_comm;
    
};

#endif /* IO_hpp */
