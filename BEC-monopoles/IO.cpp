//
//  IO.cpp
//  BEC-monopoles
//
//  Created by Adith Ramamurti on 11/9/15.
//  Copyright © 2015 Adith Ramamurti. All rights reserved.
//

#include "IO.hpp"
#include <iomanip>

void IO::set_up_outfiles(int world_rank){
    std::stringstream sstm;
    sstm << "overview_" << world_rank <<".txt";
    std::string result = sstm.str();
    f1.open(result.c_str());
    
    std::stringstream sstm2;
    sstm2 << "energydata_" << world_rank <<".csv";
    std::string result2 = sstm2.str();
    f2.open(result2.c_str());
    
    std::stringstream sstm3;
    sstm3 << "permcycledata_" << world_rank <<".csv";
    std::string result3 = sstm3.str();
    f3.open(result3.c_str());
    
    std::stringstream sstm4;
    sstm4 << "windingdata_" << world_rank <<".csv";
    std::string result4 = sstm4.str();
    f4.open(result4.c_str());
    
    f2 << "[Step No." << ",\t" << "Energy/atom" << ",\t" << "KE/atom" << ",\t" << "PE/atom" << ",\t"<< "Num. Particles]" << std::endl;
}

void IO::write_final(double energy, double energystd, int num_particles, iiVector cycles, std::vector<boost::tuple<std::string, int, int> > accept, double particles){
    
    f1 << "[Energy Statistics]" << "\n"<< std::endl;
    f1 << std::left << std::setw(25) << "Total Energy"<< std::left << std::setw(5) <<"=\t" << std::left << std::setw(15) << energy << " +/- "<< energystd<< std::endl;
    f1 << std::left << std::setw(25) << "Energy/atom" << std::left << std::setw(5) <<"=\t" << std::left << std::setw(15) << energy/num_particles <<std::endl;
    f1 << std::left << std::setw(25) << "Avg. Num. Particles" << std::left << std::setw(5) <<"=\t" << std::left << std::setw(15) << particles << "\n" << std::endl;

    int sum = 0;
    for(iiVector::iterator it = cycles.begin(); it != cycles.end(); it++){
        iVector stepcyc = *it;
        sum += std::accumulate(stepcyc.begin(), stepcyc.end(), 0);
    }

    
    f1 << "[Permutation Statistics]" <<"\n"<< std::endl;
    
    if(cycles.size() > 0){
        iVector cyclesum(cycles[0].size(),0);
        for(iiVector::iterator it = cycles.begin(); it != cycles.end(); it++){
            iVector stepcyc = *it;
            for(iVector::iterator it2 = stepcyc.begin(); it2 != stepcyc.end(); it2++){
                cyclesum[it2-stepcyc.begin()] += *it2;
            }
        }
        
        dVector cyclepercent;
        for(iVector::iterator it = cyclesum.begin(); it != cyclesum.end(); it++){
            cyclepercent.push_back(((double)*it)/sum);
        }
        for(dVector::iterator it = cyclepercent.begin(); it != cyclepercent.end(); it++){
            f1 << it-cyclepercent.begin()+1 <<":\t" << *it <<std::endl;
        }
    }
    else
        f1 << 0 <<":\t"<<0 << std::endl;
    
    f1 << std::endl;
    
    f1 << "[Acceptance Statistics]" <<"\n"<<std::endl;
    for(std::vector<boost::tuple<std::string, int, int> >::iterator it = accept.begin(); it != accept.end(); it++){
        f1 << std::left << std::setw(20) << it->get<0>() << std::left<< std::setw(6) << "--\t"<< std::left<< std::setw(8) <<"Attempts:\t" << std::left<< std::setw(10) <<it->get<2>()<<std::left<< std::setw(6) <<"--\t" << std::left<< std::setw(8) << "Accepts:\t"<< std::left<< std::setw(10) <<it->get<1>() << std::left<< std::setw(6)<<"--\t" << std::left<< std::setw(7) << "Ratio:\t" << std::left<< std::setw(10) << double(it->get<1>())/it->get<2>() << std::endl;;
    }
}


void IO::write_step_state(int step, dVector energy, iVector cycles, int num_particles, iVector wnum){
    f2 << "["<< step << "],\t" << energy[0]/num_particles << ",\t" << energy[1]/num_particles << ",\t" << energy[2]/num_particles << ",\t" << num_particles << std::endl;
    f3 << "["<< step << "],\t";
    for(iVector::iterator it = cycles.begin(); it != cycles.end(); it++){
        f3 << *it;
        if(cycles.size() - (it-cycles.begin()) != 1)
            f3 << ", ";
    }
    f3 << std::endl;
    
    f4 << "["<< step << "],\t";
    for(iVector::iterator it = wnum.begin(); it != wnum.end(); it++){
        f4 << *it;
        if(wnum.size() - (it-wnum.begin()) != 1)
            f4 << ", ";
    }
    f4 << std::endl;
}

void IO::write_parameters(boost::shared_ptr<Parameters> params){
    f1 << "[Simulation Parameters]\n";
    f1 << std::left << std::setw(25) << "Particle Type" <<std::left << std::setw(10) << ":" << std::right << std::setw(25) <<  params->get_particle_type() << "\n";
    f1 << std::left << std::setw(25) << "Num. Particles" <<std::left << std::setw(10) << "=" << std::right << std::setw(25) << params->get_num_particles() << "\n";
    f1 << std::left << std::setw(25) << "Num. Time Slices" <<std::left << std::setw(10) << "=" << std::right << std::setw(25) << params->get_num_timeslices()<< "\n";
    f1 << std::left << std::setw(25) << "Num. Dimensions" <<std::left << std::setw(10) << "=" << std::right << std::setw(25) << params->get_ndim() <<"\n";
    f1 << std::left << std::setw(25) << "Box Size" <<std::left << std::setw(10) << "=" << std::right << std::setw(25) << params->get_box_size() <<"\n";
    f1 << std::left << std::setw(25) << "tau" <<std::left << std::setw(10) << "=" << std::right << std::setw(25) << params->get_tau() << "\n";
    f1 << std::left << std::setw(25) << "lambda" <<std::left << std::setw(10) << "=" << std::right << std::setw(25) << params->get_lambda() <<"\n";
    f1 << std::left << std::setw(25) << "Temperature" <<std::left << std::setw(10) << "=" << std::right << std::setw(25) << params->get_T() << "\n" << std::endl;
}

void IO::write_equil_parameters(boost::shared_ptr<Parameters> params, double delta){
    f1 << "[Equilibration Parameters]\n";
    f1 << std::left << std::setw(25) << "Cent. of Mass delta" <<std::left << std::setw(10) << "=" << std::right << std::setw(25) << delta << "\n";
    f1 << std::left << std::setw(25) << "Chemical Potential" <<std::left << std::setw(10) << "=" << std::right << std::setw(25) << params->get_mu() <<"\n";
    f1 << std::left << std::setw(25) << "Worm Constant" <<std::left << std::setw(10) << "=" << std::right << std::setw(25) << params->get_C0() << "\n" << std::endl;
}

boost::shared_ptr<Parameters> IO::read_parameters(std::string infile){
    std::ifstream fin(infile.c_str());
    typedef boost::shared_ptr<boost::unordered_map<std::string, std::string> > sectionmap;
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
            parameters_map.insert(std::make_pair(key, sectionmap(new boost::unordered_map<std::string, std::string> )));
        }
        else if (a != ""){
            boost::unordered_map<std::string, sectionmap >::iterator it =  parameters_map.find(key);
            (*it).second->insert(std::make_pair(a, c));
        }
        
        
    }
    
    boost::shared_ptr<Parameters> params(new Parameters(parameters_map));
    return params;
}
