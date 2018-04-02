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


/* This file contains the main method for the program */
int main(int argc, char * argv[]) {

    //Initialize the MPI
    int id, numprocs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &id);
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
    int pper = numprocs; //processors per parameters file (all processors set to one parameters file by default)
    
    int start = 0; //parameters file start number (zero by default)
    std::string parameters_file;
    if(argc == 3){ //if setting pper and start
        start = std::stoi(argv[1]); //parameters file start number set by argument
        pper = std::stoi(argv[2]); //processors per parameters file set by argument
    }
    else if(argc == 2) //if setting only start
        start = std::stoi(argv[1]); //parameters file start number set by argument
    if(pper > numprocs){
        std::cerr << "Number of processors requested greater than number of processors available. Exiting..." << std::endl;
        return 1;
    }
    int world_rank = id/pper+start; //defines the group process that the individual processor belongs to (id/pper (int/int cuts off decimals) defines which parameters file and start offsets by start)
    
    std::stringstream sstm;
    sstm << "input/parameters_" << world_rank <<".cfg";
    parameters_file = sstm.str();
    
    //define the MPI sub communicator with the ids and world ranks within subprocesses
    MPI_Comm sub_comm;
    MPI_Comm_split(MPI_COMM_WORLD, id/pper, id, &sub_comm);
    int working_rank, num_workers;
    MPI_Comm_rank(sub_comm, &working_rank);
    MPI_Comm_size(sub_comm, &num_workers);
    
    
    //read in parameters for each subprocess
    IO writer(working_rank, world_rank, sub_comm);
    Parameters params;
    writer.read_parameters(parameters_file, params);
    params.num_workers = num_workers;
    
    //check for errors in parameters
    if(num_workers > params.total_slices || params.total_slices%num_workers != 0){
        std::cerr << "Number of time slices not divisible by number of worker processes. Exiting..." << std::endl;
        return 2;
    }
    if(working_rank == 0)
        std::cout << "Initializing Data Structures ... " << std::endl;
    
    //set parameters for each individual processor
    params.slices_per_process = params.total_slices/params.num_workers;
    params.my_start = (params.slices_per_process)*(working_rank);
    params.my_end = (params.slices_per_process)*(working_rank+1) - 1;
    params.set_multistep();
    
    //write out parameters used for simulation
    if(working_rank == 0)
        writer.write_parameters(params);
    
    //initialize runner for each processor and run the simulation
    Runner runner(working_rank, params, sub_comm);
    runner.run(working_rank, params, writer);
    
    //close the file writer and finalize MPI
    if(working_rank == 0){
        writer.close();
        std::cout << "Done." << std::endl;
    }
    MPI_Finalize();
    return 0;
}
