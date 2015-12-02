//
//  IO.cpp
//  BEC-monopoles
//
//  Created by Adith Ramamurti on 11/9/15.
//  Copyright © 2015 Adith Ramamurti. All rights reserved.
//

#include "IO.hpp"

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
    
    f2 << "Step No." << ", " << "Energy/atom" << ", " << "KE/atom" << ", " << "PE/atom" << std::endl;
}

void IO::write_final(double energy, double energystd, int num_particles, iiVector cycles, iVector accept){
    
    f1 << "Total Energy = " <<energy << " +/- "<< energystd<< std::endl;
    f1 << "Energy/atom= " <<energy/num_particles<< "\n" <<std::endl;
    
    int sum = 0;
    for(iiVector::iterator it = cycles.begin(); it != cycles.end(); it++){
        iVector stepcyc = *it;
        sum += std::accumulate(stepcyc.begin(), stepcyc.end(), 0);
    }
    
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
    
    f1 << "Permutation fraction" << std::endl;
    for(dVector::iterator it = cyclepercent.begin(); it != cyclepercent.end(); it++){
        f1 << it-cyclepercent.begin()+1 <<":\t" << *it <<std::endl;
    }
    
    f1 << "\nCenter of mass acceptance: "<< 1.0*accept[0]/accept.back()<< "\nBisection acceptance: "<< 1.0*accept[1]/accept.back() << "\n" << std::endl;
}


void IO::write_step_state(int step, dVector energy, iVector cycles, int num_particles, iVector wnum){
    f2 << step << ", " << energy[0]/num_particles << ", " << energy[1]/num_particles << ", " << energy[2]/num_particles << std::endl;
    for(iVector::iterator it = cycles.begin(); it != cycles.end(); it++){
        f3 << *it;
        if(cycles.size() - (it-cycles.begin()) != 1)
            f3 << ", ";
    }
    f3 << std::endl;
    
    for(iVector::iterator it = wnum.begin(); it != wnum.end(); it++){
        f4 << *it;
        if(wnum.size() - (it-wnum.begin()) != 1)
            f4 << ", ";
    }
    f4 << std::endl;
}

void IO::write_parameters(boost::shared_ptr<Parameters> params){
    f1 << "Simulation Parameters:\nN      = \t" << params->get_num_particles() << "\nNumber of Time Slices   =\t" << params->get_num_timeslices()<< "\nndim      = \t" << params->get_ndim() <<"\nBox Size      = \t" << params->get_box_size() <<"\ntau    = \t" << params->get_tau() << "\n" << "lambda =\t" << params->get_lambda() <<"\nT      = \t" << params->get_T() << "\n" << std::endl;

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
