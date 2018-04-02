//
//  moves.hpp
//
//
//  Created by Adith Ramamurti on 9/12/16.
//
//

#ifndef moves_h
#define moves_h

#include "parameters.hpp"
#include "utility.hpp"
#include "paths.hpp"

//Base class for moves (see .cpp for implementation)
class Moves{
protected:
    int num_attempts;
    int num_accepts;
    bool worm_on;
    bool worm_off;
    std::string move_name;
    MPI_Comm local_comm;
        
public:
    Moves(MPI_Comm &local);
    virtual ~Moves(){};
    virtual Moves* clone() const = 0;
    virtual double get_delta(){return 0;} //gets the delta parameter (for CoM moves)
    virtual void shift_delta(double shift){} //shifts delta
    virtual void set_delta(double shift){} //sets delta
    virtual int attempt(int &id, Parameters &params, Paths &paths, RNG &rng, Cos &cos); //attempt a move
    bool check(); //check the action of the configuration
    void accept(); //accept the move
    void reject(); //reject the move
    int get_num_accepts(){return num_accepts;}
    int get_num_attempts(){return num_attempts;}
    void reset_acceptance_counters();
    std::string get_move_name(){return move_name;}
    std::pair<bool, bool> get_move_type(){return std::pair<bool, bool>(worm_off, worm_on);}
    
    void initialize_potential_table(int &id, Parameters &params,Paths &paths, Cos &cos);
};


class Center_of_Mass : public Moves{
private:
    double delta;
    std::vector<int> ptcls;
    bool ac_re;
    std::vector<std::vector<double> > new_coordinates; //new coordinates for move
    std::vector<std::vector<double> > new_coordinates_ahead; //new coordinates for kinetic sep. table
    std::vector<int> keys_ahead; //keys for kinetic sep. table
    std::vector<std::tuple<std::pair<int, int>, std::vector<double>, double> > new_distances; //new distances for sep. table
    std::vector<std::tuple<std::pair<int, int>, double> > new_potentials; //new potential values for sep. table


public:
    Center_of_Mass(int &id, Parameters &params, MPI_Comm &local);
    ~Center_of_Mass(){};
    
    Center_of_Mass* clone() const{return new Center_of_Mass(*this);}
    
    int attempt(int &id,Parameters &params,Paths &paths, RNG &rng, Cos &cos);
    
    void shift_delta(double shift){
        delta += shift*delta;
    }
    void set_delta(double shift){
        delta = shift;
    }
    double get_delta(){
        return delta;
    }
    
    void check(int &id, Parameters &params,Paths &paths, RNG &rng, Cos &cos);
    void accept();
    void reject();
};

class Pair_Center_of_Mass : public Moves{
private:
    double delta;
    std::vector<int> ptcls;
    bool ac_re;
    std::vector<std::vector<double> > new_coordinates;
    std::vector<std::vector<double> > new_coordinates_ahead;
    std::vector<int> keys_ahead;
    std::vector<std::tuple<std::pair<int, int>, std::vector<double>, double> > new_distances;
    std::vector<std::tuple<std::pair<int, int>, double> > new_potentials;
    
    
public:
    Pair_Center_of_Mass(int &id, Parameters &params, MPI_Comm &local);
    ~Pair_Center_of_Mass(){};
    
    Pair_Center_of_Mass* clone() const{return new Pair_Center_of_Mass(*this);}
    
    int attempt(int &id,Parameters &params,Paths &paths, RNG &rng, Cos &cos);
    void shift_delta(double shift){
        delta += shift*delta;
    }
    void set_delta(double shift){
        delta = shift;
    }
    double get_delta(){
        return delta;
    }
    
    void check(int &id, Parameters &params,Paths &paths, RNG &rng, Cos &cos);
    void accept();
    void reject();
};

class Bisection : public Moves{
private:
    int multistep_dist; //distance for a bisection move
    std::vector<std::vector<int> > multisteps; //vector of possible positions in imaginary time
    std::vector<std::vector<int> > minp; //array that stores which multistep processes belong to which processor
    std::vector<int> bisection_info; //bisection info array (to be broadcast by main process)
    std::vector<int> ptcl_slice;
    bool ac_re;
    std::vector<std::vector<double> > new_coordinates;
    std::vector<std::vector<double> > new_coordinates_ahead;
    std::vector<int> keys_ahead;
    std::vector<std::tuple<std::pair<int, int>, std::vector<double>, double> > new_distances;
    std::vector<std::tuple<std::pair<int, int>, double> > new_potentials;

public:
    Bisection(int &id, Parameters &params,MPI_Comm &local);
    ~Bisection(){};
    
    Bisection* clone() const{return new Bisection(*this);}
    
    int attempt(int &id,Parameters &params,Paths &paths, RNG &rng, Cos &cos);
    void check(int &id, Parameters &params,Paths &paths,  RNG &rng, Cos &cos);
    void accept();
    void reject();
};

class Permutation_Bisection : public Moves{
private: //same as above Bisection class
    int multistep_dist;
    std::vector<std::vector<int> > multisteps;
    std::vector<std::vector<int> > minp;
    std::vector<int> bisection_info;
    std::vector<int> ptcl_slice_1;
    std::vector<int> ptcl_slice_2;
    bool ac_re;
    int choice;
    bool keep_going;
    std::vector<std::vector<double> > new_coordinates;
    std::vector<std::vector<double> > new_coordinates_ahead;
    std::vector<int> keys_ahead;
    std::vector<std::tuple<std::pair<int, int>, std::vector<double>, double> > new_distances;
    std::vector<std::tuple<std::pair<int, int>, double> > new_potentials;
    
public:
    Permutation_Bisection(int &id, Parameters &params,MPI_Comm &local);
    ~Permutation_Bisection(){};
    
    Permutation_Bisection* clone() const{return new Permutation_Bisection(*this);}
    
    int attempt(int &id,Parameters &params,Paths &paths, RNG &rng, Cos &cos);
    void check(int &id, Parameters &params,Paths &paths,  RNG &rng, Cos &cos);
    void accept();
    void reject();
};

class Open : public Moves{
private:
    bool ac_re;
    std::vector<int> open_info;
    
public:
    Open(int &id, Parameters &params, MPI_Comm &local) : Moves(local){
        move_name = "Open";
        worm_off = true;
        worm_on = false;
    };
    ~Open(){};
    
    Open* clone() const{return new Open(*this);}
    
    int attempt(int &id, Parameters &params, Paths &paths, RNG &rng, Cos &cos);
    void check(int &id, Parameters &params, Paths &paths, RNG &rng, Cos &cos);
    void accept();
    void reject();
    
};

class Close : public Moves{
private:
    bool ac_re;
    std::vector<std::vector<double> > new_coordinates;
    std::vector<std::vector<double> > new_coordinates_ahead;
    std::vector<int> keys_ahead;
    std::vector<std::vector<std::tuple<int, std::vector<double>, double> > > new_distances;
    std::vector<std::vector<std::tuple<int, double> > > new_potentials;
    
public:
    Close(int &id, Parameters &params, MPI_Comm &local) : Moves(local) {move_name = "Close";
        worm_off = false;
        worm_on = true;
    };
    ~Close(){};
    
    Close* clone() const{return new Close(*this);}
    
    int attempt(int &id, Parameters &params, Paths &paths, RNG &rng, Cos &cos);
    void check(int &id, Parameters &params, Paths &paths, RNG &rng, Cos &cos);
    void accept();
    void reject();
    
};


class Insert : public Moves{
private:
    bool ac_re;
    std::vector<int> insert_info;
    std::vector<std::vector<double> > new_coordinates;
    std::vector<std::vector<double> > new_coordinates_ahead;
    std::vector<int> keys_ahead;
    std::vector<std::vector<std::tuple<int, std::vector<double>, double> > > new_distances;
    std::vector<std::vector<std::tuple<int, double> > > new_potentials;
    
public:
    Insert(int &id, Parameters &params, MPI_Comm &local) : Moves(local) {move_name = "Insert";
        worm_off = true;
        worm_on = false;
};
    ~Insert(){};
    
    Insert* clone() const{return new Insert(*this);}
    
    int attempt(int &id, Parameters &params, Paths &paths, RNG &rng, Cos &cos);
    void check(int &id, Parameters &params, Paths &paths, RNG &rng, Cos &cos);
    void accept();
    void reject();
    
};


class Remove : public Moves{
private:
    bool ac_re;
    std::vector<int> keys_ahead;
    
public:
    Remove(int &id, Parameters &params, MPI_Comm &local) : Moves(local) {move_name = "Remove";
        worm_off = false;
        worm_on = true;
};
    ~Remove(){};
    
    Remove* clone() const{return new Remove(*this);}
    
    int attempt(int &id, Parameters &params, Paths &paths, RNG &rng, Cos &cos);
    void check(int &id, Parameters &params, Paths &paths, RNG &rng, Cos &cos);
    void accept();
    void reject();
    
};

class Advance_Tail : public Moves{
private:
    bool ac_re;
    int M;
    std::vector<std::vector<double> > new_coordinates;
    std::vector<std::vector<double> > new_coordinates_ahead;
    std::vector<int> keys_ahead;
    std::vector<std::vector<std::tuple<int, std::vector<double>, double> > > new_distances;
    std::vector<std::vector<std::tuple<int, double> > > new_potentials;
    
public:
    Advance_Tail(int &id, Parameters &params, MPI_Comm &local) : Moves(local) {move_name = "Advance Tail";
        worm_off = false;
        worm_on = true;
};
    ~Advance_Tail(){};
    
    Advance_Tail* clone() const{return new Advance_Tail(*this);}
    
    int attempt(int &id, Parameters &params, Paths &paths, RNG &rng, Cos &cos);
    void check(int &id, Parameters &params, Paths &paths, RNG &rng, Cos &cos);
    void accept();
    void reject();
    
};


class Advance_Head : public Moves{
private:
    bool ac_re;
    int M;
    std::vector<std::vector<double> > new_coordinates;
    std::vector<std::vector<double> > new_coordinates_ahead;
    std::vector<int> keys_ahead;
    std::vector<std::vector<std::tuple<int, std::vector<double>, double> > > new_distances;
    std::vector<std::vector<std::tuple<int, double> > > new_potentials;
    
public:
    Advance_Head(int &id, Parameters &params, MPI_Comm &local) : Moves(local) {move_name = "Advance Head";
        worm_off = false;
        worm_on = true;
};
    ~Advance_Head(){};
    
    Advance_Head* clone() const{return new Advance_Head(*this);}
    
    int attempt(int &id, Parameters &params, Paths &paths, RNG &rng, Cos &cos);
    void check(int &id, Parameters &params, Paths &paths, RNG &rng, Cos &cos);
    void accept();
    void reject();
    
};

class Recede_Head : public Moves{
private:
    bool ac_re;
    int M;
    std::vector<int> keys_ahead;
    
public:
    Recede_Head(int &id, Parameters &params, MPI_Comm &local) : Moves(local) {move_name = "Recede Head";
        worm_off = false;
        worm_on = true;
};
    ~Recede_Head(){};
    
    Recede_Head* clone() const{return new Recede_Head(*this);}
    
    int attempt(int &id, Parameters &params, Paths &paths, RNG &rng, Cos &cos);
    void check(int &id, Parameters &params, Paths &paths, RNG &rng, Cos &cos);
    void accept();
    void reject();
    
};

class Recede_Tail : public Moves{
private:
    bool ac_re;
    int M;
    std::vector<int> keys_ahead;
    
public:
    Recede_Tail(int &id, Parameters &params, MPI_Comm &local) : Moves(local) {move_name = "Recede Tail";
        worm_off = false;
        worm_on = true;
};
    ~Recede_Tail(){};
    
    Recede_Tail* clone() const{return new Recede_Tail(*this);}
    
    int attempt(int &id, Parameters &params, Paths &paths, RNG &rng, Cos &cos);
    void check(int &id, Parameters &params, Paths &paths, RNG &rng, Cos &cos);
    void accept();
    void reject();
    
};

class Swap_Tail : public Moves{
private:
    bool ac_re;
    double sigma_I;
    double sigma_Z;
    int choice;
    bool keep_going;
    std::vector<std::vector<double> > new_coordinates;
    std::vector<std::vector<double> > new_coordinates_ahead;
    std::vector<int> keys_ahead;
    std::vector<std::tuple<std::pair<int, int>, std::vector<double>, double> > new_distances;
    std::vector<std::tuple<std::pair<int, int>, double> > new_potentials;
    
public:
    Swap_Tail(int &id, Parameters &params, MPI_Comm &local) : Moves(local) {move_name = "Swap Tail";
        worm_off = false;
        worm_on = true;
};
    ~Swap_Tail(){};
    
    Swap_Tail* clone() const{return new Swap_Tail(*this);}
    
    int attempt(int &id, Parameters &params, Paths &paths, RNG &rng, Cos &cos);
    void check(int &id, Parameters &params, Paths &paths, RNG &rng, Cos &cos);
    void accept();
    void reject();
    
};
class Swap_Head : public Moves{
private:
    bool ac_re;
    double sigma_I;
    double sigma_Z;
    int choice;
    bool keep_going;
    std::vector<std::vector<double> > new_coordinates;
    std::vector<std::vector<double> > new_coordinates_ahead;
    std::vector<int> keys_ahead;
    std::vector<std::tuple<std::pair<int, int>, std::vector<double>, double> > new_distances;
    std::vector<std::tuple<std::pair<int, int>, double> > new_potentials;
    
public:
    Swap_Head(int &id, Parameters &params, MPI_Comm &local) : Moves(local) {move_name = "Swap Head";        worm_off = false;
        worm_on = true;
};
    ~Swap_Head(){};
    
    Swap_Head* clone() const{return new Swap_Head(*this);}
    
    int attempt(int &id, Parameters &params, Paths &paths, RNG &rng, Cos &cos);
    void check(int &id, Parameters &params, Paths &paths, RNG &rng, Cos &cos);
    void accept();
    void reject();
};

#endif /* moves_h */
