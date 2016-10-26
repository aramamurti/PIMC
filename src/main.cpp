//
//  main.cpp
//  PIMC
//
//  Created by Adith Ramamurti on 10/13/16.
//  Copyright © 2016 Adith Ramamurti. All rights reserved.
//

#include <iostream>
#include <mpi.h>
#include <sstream>
#include "IO.hpp"
#include "parameters.hpp"
#include "runner.hpp"

int main(int argc, char * argv[]) {
    int id, numprocs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &id);
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);

    int pper = numprocs;
    int start = 0;

    std::string parameters_file;
    if(argc == 3){
        pper = std::stoi(argv[1]);
        start = std::stoi(argv[2]);
    }
    else if(argc == 2)
        pper = std::stoi(argv[1]);

    int world_rank = id/pper+start;
    std::stringstream sstm;
    sstm << "input/parameters_" << world_rank <<".cfg";
    parameters_file = sstm.str();
    
    MPI_Comm sub_comm;
    MPI_Comm_split(MPI_COMM_WORLD, id/pper, id, &sub_comm);

    int working_rank, num_workers;
    MPI_Comm_rank(sub_comm, &working_rank);
    MPI_Comm_size(sub_comm, &num_workers);
    IO writer(working_rank, world_rank, sub_comm);
    Parameters params;
    writer.read_parameters(parameters_file, params);
    params.num_workers = num_workers;
    params.slices_per_process = params.total_slices/params.num_workers;
    params.my_start = (params.slices_per_process)*(working_rank);
    params.my_end = (params.slices_per_process)*(working_rank+1) - 1;
    if(working_rank == 0)
        writer.write_parameters(params);
    Runner runner(working_rank, params, sub_comm);
    runner.run(working_rank, params, writer);
    if(working_rank == 0)
        writer.close();
    MPI_Finalize();
    return 0;
}
