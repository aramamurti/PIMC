//
//  IO.cpp
//  BEC-monopoles
//
//  Created by Adith Ramamurti on 11/9/15.
//  Copyright © 2015 Adith Ramamurti. All rights reserved.
//

#include "IO.hpp"
#include <iomanip>

void IO::set_up_outfiles(int& id){
    std::stringstream sstm;
    sstm << "output/overview_" << id <<".txt";
    std::string result = sstm.str();
    f1.open(result.c_str());
    
    sstm.str(std::string());
    sstm.clear();
    sstm << "output/parameters_" << id <<".csv";
    result = sstm.str();
    f5.open(result.c_str());
    
    sstm.str(std::string());
    sstm.clear();
    sstm << "output/energy_data_" << id <<".csv";
    result = sstm.str();
    f2.open(result.c_str());
    
    sstm.str(std::string());
    sstm.clear();
    sstm << "output/permutation_data_" << id <<".csv";
    result = sstm.str();
    f3.open(result.c_str());
    
    sstm.str(std::string());
    sstm.clear();
    sstm << "output/winding_data_" << id <<".csv";
    result = sstm.str();
    f4.open(result.c_str());
    
    f2 << "[Step No." << ",\t" << "TE/atom (th.)" << ",\t" << "KE/atom (th.)" << ",\t" << "PE/atom (th.)" << ",\t"<< "TE/atom (v.)" << ",\t" << "KE/atom (v.)" << ",\t" << "PE/atom (v.)" << ",\t"<<"Num. Particles]" << std::endl;

    sstm.str(std::string());
    sstm.clear();
    sstm << "output/path_dump_" << id <<".txt";
    result = sstm.str();
    pathf.open(result.c_str());
    pathf<<"{";
    pathf.flush();
    
    sstm.str(std::string());
    sstm.clear();
    sstm << "output/connections_" << id <<".txt";
    result = sstm.str();
    pathf2.open(result.c_str());
    pathf2 << "{";
    pathf2.flush();
    
    sstm.str(std::string());
    sstm.clear();
    sstm << "output/charges_" << id <<".txt";
    result = sstm.str();
    pathf3.open(result.c_str());
    pathf3 << "{";
    pathf3.flush();

}

void IO::write_final(int counter, double particles, double energy_th, double energy_th_std, double energy_v, double energy_v_std, std::vector<std::vector<int> >& cycles, boost::ptr_vector<Moves> &moves){
    f1 << std::left << std::setw(25) << "Num. Samples" << std::left << std::setw(5) <<"=\t" << std::left << std::setw(15) << counter << "\n" << std::endl;
    f1 << std::endl;
    f1 << "[Energy Statistics]" << "\n"<< std::endl;
    f1 << std::left << std::setw(25) << "Avg. Num. Particles" << std::left << std::setw(5) <<"=\t" << std::left << std::setw(15) << particles << "\n" << std::endl;
    f1 << std::left << std::setw(25) << "Total Energy (Th.)"<< std::left << std::setw(5) <<"=\t" << std::left << std::setw(15) << energy_th << " +/- "<< energy_th_std<< std::endl;
    f1 << std::left << std::setw(25) << "Total Energy (V.)"<< std::left << std::setw(5) <<"=\t" << std::left << std::setw(15) << energy_v << " +/- "<< energy_v_std<< std::endl;
    f1 << std::left << std::setw(25) << "Energy/Particle (Th.)"<< std::left << std::setw(5) <<"=\t" << std::left << std::setw(15) << energy_th/particles << " +/- "<< energy_th_std/particles << std::endl;
    f1 << std::left << std::setw(25) << "Energy/Particle (V.)"<< std::left << std::setw(5) <<"=\t" << std::left << std::setw(15) << energy_v/particles << " +/- "<< energy_v_std/particles << std::endl;
    int sum = 0;
    for(std::vector<std::vector<int> >::iterator it = cycles.begin(); it != cycles.end(); ++it){
        std::vector<int> stepcyc = *it;
        sum += std::accumulate(stepcyc.begin(), stepcyc.end(), 0);
    }
    f1 << "[Permutation Statistics]" <<"\n"<< std::endl;
    if(cycles.size() > 0){
        std::vector<int> cyclesum(ceil(particles),0);
        for(std::vector<std::vector<int> >::iterator it = cycles.begin(); it != cycles.end(); ++it){
            std::vector<int> stepcyc = *it;
            if(stepcyc.size() > cyclesum.size())
                cyclesum.resize(stepcyc.size());
            for(std::vector<int>::iterator it2 = stepcyc.begin(); it2 != stepcyc.end(); ++it2){
                cyclesum[it2-stepcyc.begin()] += *it2;
            }
        }
        std::vector<double> cyclepercent;
        for(std::vector<int>::iterator it = cyclesum.begin(); it != cyclesum.end(); ++it){
            cyclepercent.push_back(((double)*it)/sum);
        }
        for(std::vector<double>::iterator it = cyclepercent.begin(); it != cyclepercent.end(); ++it){
            f1 << it-cyclepercent.begin()+1 <<":\t" << *it <<std::endl;
        }
    }
    else
        f1 << 0 <<":\t"<<0 << std::endl;
    f1 << std::endl;
    f1 << "[Acceptance Statistics]" <<"\n"<<std::endl;
    for(auto &it : moves){
        f1 << std::left << std::setw(20) << it.get_move_name() << std::left<< std::setw(6) << "--\t"<< std::left<< std::setw(8) <<"Attempts:\t" << std::left<< std::setw(10) <<it.get_num_attempts()<<std::left<< std::setw(6) <<"--\t" << std::left<< std::setw(8) << "Accepts:\t"<< std::left<< std::setw(10) << it.get_num_accepts() << std::left<< std::setw(6)<<"--\t" << std::left<< std::setw(7) << "Ratio:\t" << std::left<< std::setw(10) << double(it.get_num_accepts())/it.get_num_attempts() << std::endl;;
    }
}


void IO::write_step_state(int step, int num_particles, std::vector<double>& energy, std::vector<int>& cycles, std::vector<int>& wnum){
    f2 << "["<< step << "],\t" << energy[0]/num_particles << ",\t" << energy[1]/num_particles << ",\t" << energy[2]/num_particles << ",\t" << energy[3]/num_particles << ",\t" << energy[4]/num_particles << ",\t" << energy[5]/num_particles << ",\t"<< num_particles << std::endl;
    f3 << "["<< step << "],\t";
    for(std::vector<int>::iterator it = cycles.begin(); it != cycles.end(); ++it){
        f3 << *it;
        if(cycles.size() - (it-cycles.begin()) != 1)
            f3 << ", ";
    }
    f3 << std::endl;
    
    f4 << "["<< step << "],\t";
    for(std::vector<int>::iterator it = wnum.begin(); it != wnum.end(); ++it){
        f4 << *it;
        if(wnum.size() - (it-wnum.begin()) != 1)
            f4 << ", ";
    }
    f4 << std::endl;
}

void IO::write_parameters(Parameters& params){
    f1 << "[Simulation Parameters]\n";
    f1 << std::left << std::setw(25) << "Particle Type" <<std::left << std::setw(10) << ":" << std::right << std::setw(25) <<  params.particle_type << "\n";
    f1 << std::left << std::setw(25) << "Particles" <<std::left << std::setw(10) << "=" << std::right << std::setw(25) << params.particles << "\n";
    f1 << std::left << std::setw(25) << "Time Slices" <<std::left << std::setw(10) << "=" << std::right << std::setw(25) << params.total_slices<< "\n";
    f1 << std::left << std::setw(25) << "Dimensions" <<std::left << std::setw(10) << "=" << std::right << std::setw(25) << params.dimensions <<"\n";
    f1 << std::left << std::setw(25) << "Box Size" <<std::left << std::setw(10) << "=" << std::right << std::setw(25) << params.box_size <<"\n";
    f1 << std::left << std::setw(25) << "tau" <<std::left << std::setw(10) << "=" << std::right << std::setw(25) << params.tau << "\n";
    f1 << std::left << std::setw(25) << "lambda" <<std::left << std::setw(10) << "=" << std::right << std::setw(25) << params.lambda <<"\n";
    f1 << std::left << std::setw(25) << "Temperature" <<std::left << std::setw(10) << "=" << std::right << std::setw(25) << params.temperature << "\n";
    f1 << std::left << std::setw(25) << "Coupling" <<std::left << std::setw(10) << "=" << std::right << std::setw(25) << params.coupling << "\n" << std::endl;

    f5 << "Particle Type, " << params.particle_type << std::endl;
    f5 << "Particles, " << params.particles << std::endl;
    f5 << "Dimensions, " << params.dimensions << std::endl;
    f5 << "Temperature, " << params.temperature << std::endl;
    f5 << "Time Slices, " << params.total_slices << std::endl;
    f5 << "Box Size, " << params.box_size << std::endl;
    f5 << "Tau, " << params.tau << std::endl;
    f5 << "Lambda, " << params.lambda << std::endl;
    f5 << "Coupling, " << params.coupling << std::endl;
}

void IO::write_equil_parameters(Parameters& params, double delta, double delta_pair){
    f1 << "[Equilibration Parameters]\n";
    f1 << std::left << std::setw(25) << "Cent. of Mass delta" <<std::left << std::setw(10) << "=" << std::right << std::setw(25) << delta << "\n";
    if(delta_pair != 0)
        f1 << std::left << std::setw(25) << "Pair Cent. of Mass delta" <<std::left << std::setw(10) << "=" << std::right << std::setw(25) << delta_pair << "\n";
    f1 << std::left << std::setw(25) << "C0" <<std::left << std::setw(10) << "=" << std::right << std::setw(25) << params.C0 << "\n";
    f1 << std::left << std::setw(25) << "µ" <<std::left << std::setw(10) << "=" << std::right << std::setw(25) << params.mu << "\n";
    f1 << std::endl;
}

void IO::read_parameters(std::string infile, Parameters& params){
    std::ifstream fin(infile.c_str());
    typedef boost::unordered_map<std::string, std::string> sectionmap;
    boost::unordered_map<std::string, sectionmap > parameters_map;
    
    std::string key;
    std::string line;
    while (std::getline(fin, line))
    {
        std::istringstream sstr1(line);
        std::string a, b, c;
        sstr1 >> a >> b >> c;
        if(a[0] == '['){
            std::stringstream sstr2;
            sstr2 << a.substr(1) << b;
            key = sstr2.str().substr(0,sstr2.str().length()-1);
            sectionmap section;
            parameters_map.insert(std::make_pair(key, section));
        }
        else if (a != ""){
            boost::unordered_map<std::string, sectionmap >::iterator it =  parameters_map.find(key);
            (*it).second.insert(std::make_pair(a, c));
        }
    }
    params.set_parameters(parameters_map);
}

void IO::path_dump(int& id, int count, Paths& paths, Parameters &params){
    std::vector<std::vector<std::vector<double> > > all_locations;
    if(id == 0){
        std::vector<std::vector<double> > slice;
        for(int i = 0; i < params.slices_per_process; ++i){
            for(int k = 0; k < params.particles; ++k){
                slice.push_back(paths.get_coordinate(i,k));
            }
            all_locations.push_back(slice);
            slice.clear();
        }
        std::vector<double> loc(params.dimensions);
        for(int j = 1; j < params.num_workers; ++j){
            for(int i = 0; i < params.slices_per_process; ++i){
                for(int k = 0; k < params.particles; ++k){
                    MPI_Recv(&loc[0], params.dimensions, MPI_DOUBLE, j, i*params.particles+k, local_comm, MPI_STATUS_IGNORE);
                    slice.push_back(loc);
                }
                all_locations.push_back(slice);
                slice.clear();
            }
        }
    }
    else{
        for(int i = 0; i < params.slices_per_process; ++i)
            for(int ptcl = 0; ptcl < params.particles; ++ptcl)
                MPI_Send(&paths.get_coordinate(i,ptcl)[0], params.dimensions, MPI_DOUBLE, 0, i*params.particles+ptcl, local_comm);
    }
    if(id == 0){
        pathf << std::fixed << std::showpoint;
        pathf << std::setprecision(6);
        pathf<<"{";
        for(int slice = 0; slice < params.total_slices; ++slice){
            pathf << "{";
            for(int ptcl = 0; ptcl < params.particles; ++ptcl){
                pathf << "{";
                for(int dim = 0; dim < params.dimensions; ++dim){
                    pathf << all_locations[slice][ptcl][dim];
                    if(dim != params.dimensions -1)
                        pathf << ", ";
                }
                if(ptcl != params.particles-1)
                    pathf << "},";
                else
                    pathf << "}";
            }
            pathf << "},";
        }
        pathf << "{";
        for(int ptcl = 0; ptcl < params.particles; ++ptcl){
            pathf << "{";
            for(int dim = 0; dim < params.dimensions; ++dim){
                pathf << all_locations[0][paths.forward_connects[ptcl]][dim];
                if(dim != params.dimensions -1)
                    pathf << ", ";
            }
            if(ptcl != params.particles-1)
                pathf << "},";
            else
                pathf << "}";
        }
        pathf << "}},"  ;
        pathf.flush();
        pathf2 << "{";
        for(int ptcl = 0; ptcl < params.particles; ++ptcl){
            pathf2 << paths.forward_connects[ptcl];
            if(ptcl != params.particles-1)
                pathf2 << ",";
            else
                pathf2 << "},";
        }
        pathf2.flush();
        pathf3 << "{";
        for(int ptcl = 0; ptcl < params.particles; ++ptcl){
            pathf3 << paths.charge[ptcl];
            if(ptcl != params.particles-1)
                pathf3 << ",";
            else
                pathf3 << "},";
        }
        pathf3.flush();
    }
}
void IO::write_acceptances(int counter, boost::ptr_vector<Moves> &moves){
    f1 << std::left << std::setw(25) << "Num. Samples" << std::left << std::setw(5) <<"=\t" << std::left << std::setw(15) << counter << "\n" << std::endl;
    f1 << std::endl;
    f1 << "[Acceptance Statistics]" <<"\n"<<std::endl;
    for(auto &it : moves){
        f1 << std::left << std::setw(20) << it.get_move_name() << std::left<< std::setw(6) << "--\t"<< std::left<< std::setw(8) <<"Attempts:\t" << std::left<< std::setw(10) <<it.get_num_attempts()<<std::left<< std::setw(6) <<"--\t" << std::left<< std::setw(8) << "Accepts:\t"<< std::left<< std::setw(10) << it.get_num_accepts() << std::left<< std::setw(6)<<"--\t" << std::left<< std::setw(7) << "Ratio:\t" << std::left<< std::setw(10) << double(it.get_num_accepts())/it.get_num_attempts() << std::endl;;
    }
}
