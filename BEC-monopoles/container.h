//
//  beadContainer.h
//  BEC-monopoles
//
//  Created by Adith Ramamurti on 6/8/15.
//  Copyright (c) 2015 Adith Ramamurti. All rights reserved.
//

#ifndef BEC_monopoles_beadContainer_h
#define BEC_monopoles_beadContainer_h
#include "lookuptable.hpp"
#include "utility.h"
#include <iterator>

template<class T>

class PathList {
private:
    typedef std::vector<int> iVector;
    typedef std::vector<double> dVector;
    typedef std::vector<std::vector<double> > ddVector;
    typedef std::vector<std::vector<int> > iiVector;
    
    
    class Node {
    public:
        typedef boost::shared_ptr<Node> node_ptr;
        
        int row_number = 0;
        int column_number = 0;
        
        int old_row_number;
        size_t key;
        
        T data;
        T old_data;
        
        int charge;
        
        node_ptr left_node;
        node_ptr right_node;
        node_ptr old_right_node;
        node_ptr old_left_node;
        
        Node(T e) {
            this->data = e;
            left_node = node_ptr();
            right_node = node_ptr();
        }
        ~Node() {
            left_node.reset();
            right_node.reset();
        }
    };
    
    typedef boost::shared_ptr<Node> node_ptr;
    
    bool circular = false;
    
    size_t push_counter = 0;
    iVector size;
    iVector old_size;
    
    std::vector<node_ptr> head;
    std::vector<node_ptr> tail;
    
    iVector permute_start_orders;
    iVector permute_end_order;
    iVector repermute_start_order;
    iVector repermute_end_order;
    
    iVector pso;
    iVector peo;
    iVector rso;
    iVector reo;
    
    int permute_position = 0;
    
    std::vector<std::vector<node_ptr> > list_map;
    std::vector<std::vector<node_ptr> > old_list_map;
    
    int num_charges;
    
    boost::unordered_map<size_t, node_ptr> key_map;
    boost::unordered_map<size_t, node_ptr> old_key_map;
    
    boost::shared_ptr<Separation_Table> sep;
    boost::shared_ptr<Permutation_Separation_Table> perm_sep;
    
    boost::shared_ptr<Utility> util;
    
    double boxsize;
    
public:
    
    PathList(int ndim, int multistep_dist, double boxsize = -1, int num_charges = 0) {
        size.resize(0);
        sep = boost::shared_ptr<Separation_Table>(new Separation_Table(boxsize));
        perm_sep = boost::shared_ptr<Permutation_Separation_Table>(new Permutation_Separation_Table(multistep_dist,boxsize));
        util = boost::shared_ptr<Utility>(new Utility(0));
        this->num_charges = num_charges;
        this->boxsize = boxsize;
    }
    
    PathList(std::string infile)
    {
        size.resize(0);
        
        std::ifstream in(infile.c_str());
        std::string s;
        size_t numbreaks = 0;
        size_t numline = 0;
        std::vector<ddVector > dataset;
        while(!in.eof()){
            std::getline(in,s);
            if(s.length()!= 0){
                std::string delimiter = "), ";
                size_t pos = 0;
                std::string token;
                
                while ((pos = s.find(delimiter)) != std::string::npos) {
                    token = s.substr(s.find("(")+1, pos-1);
                    s.erase(0, pos+delimiter.length());
                    if(numbreaks == 0){
                        size_t n = std::count(token.begin(),token.end(),',')+1;
                        dVector data(n);
                        std::vector<std::string> ent;
                        boost::split(ent, token, boost::is_any_of(", "), boost::token_compress_on);
                        for(std::vector<std::string>::iterator it = ent.begin(); it!=ent.end(); it++){
                            data[it-ent.begin()]= std::atof((*it).c_str());
                        }
                        if(dataset.size() < numline+1)
                            dataset.resize(numline+1);
                        dataset[numline].push_back(data);
                    }
                    if(numbreaks == 1){
                        std::vector<std::string> ent;
                        boost::split(ent, token, boost::is_any_of(", "), boost::token_compress_on);
                        ddVector dset2 = dataset[numline];
                        dVector nodedat= dset2[std::atoi(ent[1].c_str())];
                        push_back(nodedat, numline, std::atoi(ent[0].c_str()),std::atoi(ent[1].c_str()));
                    }
                }
                token = s.substr(1,s.length()-2);
                if(numbreaks == 0){
                    size_t n = std::count(token.begin(),token.end(),',')+1;
                    dVector data(n);
                    std::vector<std::string> ent;
                    boost::split(ent, token, boost::is_any_of(", "), boost::token_compress_on);
                    for(std::vector<std::string>::iterator it = ent.begin(); it!=ent.end(); it++){
                        data[it-ent.begin()]= std::atof((*it).c_str());
                    }
                    if(dataset.size() < numline+1)
                        dataset.resize(numline+1);
                    dataset[numline].push_back(data);
                }
                if(numbreaks == 1){
                    std::vector<std::string> ent;
                    boost::split(ent, token, boost::is_any_of(", "), boost::token_compress_on);
                    ddVector dset2 = dataset[numline];
                    dVector nodedat= dset2[std::atoi(ent[1].c_str())];
                    push_back(nodedat, numline, std::atoi(ent[0].c_str()),std::atoi(ent[1].c_str()));
                }
                
                if(numbreaks == 2){
                    std::vector<std::string> ent;
                    boost::split(ent, s, boost::is_any_of("("), boost::token_compress_on);
                    std::string ss = ent.back();
                    ss = ss.substr(0,s.length()-2);
                    boost::split(ent, ss, boost::is_any_of(", "), boost::token_compress_on);
                    tail[numline]->right_node = head[std::atoi(ent[0].c_str())];
                    tail[numline]->right_node->left_node = tail[numline];
                }
                numline++;
            }
            else{
                numbreaks++;
                numline=0;
            }
            
        }
        circular = true;
    }
    
    ~PathList() {
        make_circular(false);
        for(typename std::vector< std::vector<node_ptr>>::iterator it = list_map.begin(); it != list_map.end(); it++){
            node_ptr temp = (*it)[0];
            while (temp != NULL) {
                node_ptr next(temp->right_node);
                temp->right_node.reset();
                temp->left_node.reset();
                temp.reset();
                temp = next;
            }
        }
    }
    
    void push_back(T t, int index = 0, size_t key = -1, int row_number = -1, int column_number = -1) {
        node_ptr node(new Node(t));
        if(row_number == -1)
            node->row_number = index;
        else{
            node->row_number = row_number;
            index = row_number;
        }
        
        if(size.size() <= index){
            size.resize(index+1);
            size[index] = 0;
        }
        ++size[index];
        
        list_map.resize(size.size());
        list_map[index].resize(size[index]);
        
        if(column_number == -1)
            node->column_number = size[index]-1;
        else
            node->column_number = column_number;
        
        if(key == -1){
            set_key_add_sep(node);
        }
        else{
            node->key = key;
        }
        
        update_key_map(node);
        
        
        if (size[index] == 1) {
            list_map[index][size[index]-1] = node;
            return;
        }
        
        node->left_node = list_map[index][size[index]-2];
        node->left_node->right_node = node;
        
        list_map[index][size[index]-1] = node;
        
        node.reset();
    }
    
    void set_key_add_sep(node_ptr node){
        node->key = push_counter;
        push_counter++;
        sep->add_bead(node->column_number, node->key, node->data);
        perm_sep->add_bead(node->column_number, node->key, node->data);
    }
    
    void set_charge(int ptcl, int charge){
        for(int slice = 0; slice < size[0]; slice++)
            list_map[ptcl][slice]->charge = charge;
    }
    
    int get_charge(int ptcl){
        return list_map[ptcl][0]->charge;
    }
    
    void make_circular(bool circ = true){
        
        if(circ){
            for(typename std::vector< std::vector<node_ptr> >::iterator it = list_map.begin(); it != list_map.end(); it++){
                std::vector<node_ptr> list = (*it);
                list[0]->left_node = list[list_map[it-list_map.begin()].size()-1];
                list[0]->left_node->right_node = list[0];
            }
            circular = true;
        }
        else{
            for(typename std::vector< std::vector<node_ptr>>::iterator it = list_map.begin(); it != list_map.end(); it++){
                node_ptr ln = (*it)[0];
                ln->left_node->right_node.reset();
                ln->left_node.reset();
            }
            circular = false;
        }
    }
    
    void update_key_map(node_ptr node){
        size_t key = node->key;
        auto it = key_map.find(key);
        if(it != key_map.end()){
            key_map.erase(it);
        }
        key_map.insert(std::pair<size_t, node_ptr>(key, node));
        
    }
    
    void remove_from_key_map(node_ptr node){
        size_t key = node->key;
        auto it = key_map.find(key);
        if(it != key_map.end())
            key_map.erase(it);
    }
    
    void generate_perm_seps(){
        for(int i = 0 ;i < list_map.size(); i++)
            for(int j = 0; j <list_map[i].size(); j++)
                perm_sep->update_bead(list_map[i][j]->column_number,list_map[i][j]->key);
    }
    
    void generate_bead_seps(){
        for(int i = 0 ;i < list_map.size(); i++)
            for(int j = 0; j <list_map[i].size(); j++)
                sep->update_bead_separations(list_map[i][j]->column_number,list_map[i][j]->key);
    }
    
    void shift_all(int row, T shift, double box_size){
        
        int times_through = 0;
        
        int cur_row = row;
        while(!(cur_row == row && times_through != 0)){
            for(int i = 0; i < size[cur_row]; i++){
                std::transform(list_map[cur_row][i]->data.begin(), list_map[cur_row][i]->data.end(), shift.begin(),
                               list_map[cur_row][i]->data.begin(), std::plus<double>());
                util->put_in_box(list_map[cur_row][i]->data, box_size);
                sep->update_bead_separations(i, list_map[cur_row][i]->key);
                perm_sep->update_bead(i, list_map[cur_row][i]->key);

            }
            cur_row = list_map[cur_row][size[cur_row]-1]->right_node->row_number;
            times_through++;
        }
    }
    
    iVector get_changed_ptcls(int row){
        iVector shifted_rows(1,row);
        int start_row = row;
        while(list_map[row][size[row]-1]->right_node->row_number != row){
            row = list_map[row][size[row]-1]->right_node->row_number;
            if(row == start_row)
                break;
            shifted_rows.push_back(row);
        }
        
        return shifted_rows;
    }
    
    void set_old_data(){
        for(int i = 0; i < size.size(); i++){
            for(int j = 0; j<size[i];j++){
                list_map[i][j]->old_data = list_map[i][j]->data;
                list_map[i][j]->old_row_number = list_map[i][j]->row_number;
                list_map[i][j]->old_right_node = list_map[i][j]->right_node;
                list_map[i][j]->old_left_node = list_map[i][j]->left_node;
            }
        }
    }
    
    void confirm_sep_update(){
        sep->confirm_update();
        perm_sep->confirm_update();
        sep->set_update(true);
        perm_sep->set_update(true);
    }
    
    void clear_old_nodes(){
        for(int i = 0; i < size.size(); i++){
            for(int j = 0; j<size[i];j++){
                list_map[i][j]->old_right_node.reset();
                list_map[i][j]->old_left_node.reset();
            }
        }
    }
    
    void revert_old_data(){
        for(int i = 0; i < size.size(); i++){
            for(int j = 0; j<size[i];j++){
                list_map[i][j]->data = list_map[i][j]->old_data;
                list_map[i][j]->row_number = list_map[i][j]->old_row_number;
                list_map[i][j]->right_node = list_map[i][j]->old_right_node;
                list_map[i][j]->left_node = list_map[i][j]->old_left_node;
                list_map[i][j]->old_right_node.reset();
                list_map[i][j]->old_left_node.reset();
            }
        }
    }
    
    void set_permutation(iVector& i, iVector& j, int pos, int dist = 0){
        permute_position = (pos+dist);
        permute_start_orders.resize(i.size());
        permute_end_order.resize(j.size());
        
        std::copy(i.begin(), i.end(), permute_start_orders.begin());
        std::copy(j.begin(), j.end(), permute_end_order.begin());
        
        if(permute_position > size[0]){
            circular_perm(i,j, repermute_start_order, repermute_end_order);
        }
        else{
            repermute_start_order = i;
            repermute_end_order = j;
        }
    }
    
    void permute(bool reverse = false){
        iVector i;
        iVector j;
        i.resize(permute_start_orders.size());
        j.resize(permute_end_order.size());
        std::copy(permute_start_orders.begin(), permute_start_orders.end(), i.begin());
        std::copy(permute_end_order.begin(), permute_end_order.end(), j.begin());
        
        
        if(!reverse){
            
            if(!i.empty()) {
                
                old_list_map = list_map;
                
                int pos = permute_position;
                if(permute_position > size[0]){
                    i = repermute_start_order;
                    j = repermute_end_order;
                    pos = pos%size[0];
                }
                
                std::vector<node_ptr> temps;
                
                for(int ptc = 0; ptc < size.size(); ptc ++)
                    temps.push_back(list_map[ptc][pos-1]);
                
                std::vector<node_ptr> temps2(j.size());
                for(typename std::vector<node_ptr>::iterator it = temps2.begin(); it != temps2.end(); it++)
                    *it = temps[j[it-temps2.begin()]]->right_node;
                
                for(typename iVector::iterator it = i.begin(); it != i.end(); it++){
                    if(permute_position != 0){
                        temps[*it]->old_right_node = temps[*it]->right_node;
                        temps[*it]->right_node = temps2[it-i.begin()];
                        temps[*it]->right_node->left_node = temps[*it];
                    }
                }
                
                for(int curpos = pos; curpos < size[0]; curpos++){
                    for(iVector::iterator it = i.begin(); it != i.end();it++){
                        list_map[*it][curpos] = old_list_map[j[it-i.begin()]][curpos];
                    }
                }
                
                reset_indices();
                
            }
        }
        else{
            if(!i.empty()) {
                int pos = permute_position;
                i = repermute_start_order;
                if(permute_position > size[0]){
                    pos = pos%size[0];
                }
                
                std::vector<node_ptr> temps;
                
                for(int ptc = 0; ptc < size.size(); ptc ++)
                    temps.push_back(list_map[ptc][pos-1]);
                
                for(typename iVector::iterator it = i.begin(); it != i.end(); it++){
                    if(permute_position != 0){
                        temps[*it]->right_node = temps[*it]->old_right_node;
                        temps[*it]->right_node->left_node = temps[*it];
                    }
                }
                
                list_map = old_list_map;
                
                reset_indices();
                
            }
        }
    }
    
        
    void circular_perm(iVector& i, iVector& j, iVector& ip, iVector& jp){
        ip.clear();
        jp.clear();
        ip.reserve(i.size());
        jp.reserve(j.size());

        for(iVector::iterator it = i.begin(); it != i.end(); it++){
            ip.push_back(list_map[*it][size[*it]-1]->right_node->row_number);
            jp.push_back(list_map[j[it-i.begin()]][size[*it]-1]->right_node->row_number);
            
        }
    }

    
    void get_perm_seps(iVector& i, iVector& j, iVector& ip, iVector& jp, std::vector<std::pair<double, double> >& pair_dists, int pos, int dist = 0){
        
        if(pos+dist > size[0]){
            circular_perm(i,j, ip, jp);
        }
        else{
            ip = i;
            jp = j;
        }
        
        int end = (pos+dist)%size[0];
        for(iVector::iterator it = i.begin(); it!=i.end(); it++){
            if(*it != j[it-i.begin()]){
                std::pair<size_t,size_t> pair_1(list_map[*it][pos]->key, list_map[jp[it-i.begin()]][end]->key);
                std::pair<size_t,size_t> pair_2(list_map[*it][pos]->key, list_map[ip[it-i.begin()]][end]->key);
                pair_dists.push_back(std::pair<double, double>(perm_sep->get_perm_separation(pair_1),perm_sep->get_perm_separation(pair_2)));
            }
        }
        
    }
    
    T get_bead_data(int row, int slice){
        
        if(slice < size[row]){
            return list_map[row][slice]->data;
        }
        else if(circular){
            slice = slice%size[row];
            row = list_map[row][size[row]-1]->right_node->row_number;
            return list_map[row][slice]->data;
        }
        else
            return T(0);
    }
    
    void set_bead_data(int row, int slice, double box_size, T data, T old_data = T(0)){
        if(slice < size[row]){
        }
        else if(circular){
            slice = slice%size[row];
            row = list_map[row][size[row]-1]->right_node->row_number;
        }
        util->put_in_box(data, box_size);
        list_map[row][slice]->data = data;
        
        sep->update_bead_separations(slice, list_map[row][slice]->key);
        perm_sep->update_bead(slice, list_map[row][slice]->key);
        
        if(old_data.size() != 0)
            list_map[row][slice]->old_data = old_data;
    }
    
    void get_pair_same_path(int row, int start, int dist, T& bead1, T& bead2){
        
        int row2 = row;
        int end = start+dist;
        
        if(start >= size[row] && end >= size[row] && circular){
            start = start%size[row];
            row = list_map[row][size[row]-1]->right_node->row_number;
        }
        
        if(end >= size[row] && circular){
            end = end%size[row];
            row2 = list_map[row2][size[row]-1]->right_node->row_number;
        }
        
        bead1 = list_map[row][start]->data;
        bead2 = list_map[row2][end]->data;
    }
    
    iVector get_pair_rows(int row, int start, int dist){
        
        int row2 = row;
        int end = start+dist;
        
        if(start >= size[row] && end >= size[row] && circular){
            start = start%size[row];
            row = list_map[row][size[row]-1]->right_node->row_number;
        }
        
        if(end >= size[row] && circular){
            end = end%size[row];
            row2 = list_map[row2][size[row]-1]->right_node->row_number;
        }
        
        iVector ret(0);
        ret.push_back(row);
        ret.push_back(row2);
        return ret;
    }
    
        
    const T& get_path_separation(int row1, int row2, int slice){
        return sep->get_bead_separation(std::pair<size_t,size_t>(list_map[row1][slice]->key,list_map[row2][slice]->key));
    }
    
    
    iVector get_cycles(){
        iVector cyclenum(size.size(),0);
        for(int row = 0; row < size.size(); row++){
            node_ptr temp = list_map[row][0]->right_node;
            int tempPos = 1;
            while(temp != list_map[row][0]){
                temp = temp->right_node;
                ++tempPos;
            }
            int cycNum = tempPos/size[row];
            
            cyclenum[cycNum-1]++;
        }
        
        return cyclenum;
    }
    
    
    iVector get_rep_swap_order(){
        return rso;
    }
    
    void set_prev_perm(){
        pso = permute_start_orders;
        peo = permute_end_order;
        rso = repermute_start_order;
        reo = repermute_end_order;
    }
    
    void reset_permute(){
        permute_start_orders = pso;
        permute_end_order = peo;
        repermute_start_order = rso;
        repermute_end_order = reo;
    }
    
    int get_num_particles(){
        return size.size();
    }
    
    
    void reset_indices(){
        for(int row = 0; row < size.size(); row++){
            for(int slice = 0; slice < size[row]; slice++){
                list_map[row][slice]->row_number = row;
            }
        }
    }
    
    bool check_indices(){
        for(int row = 0; row < size.size(); row++){
            for(int slice = 0; slice < size[row]; slice++){
                if(list_map[row][slice]->row_number != row)
                    return false;
            }
        }
        return true;
    }
    
    bool check_data(){
        for(int row = 0; row < size.size(); row++){
            for(int slice = 0; slice < size[row]; slice++){
                if(isnan(list_map[row][slice]->data[0])||isnan(list_map[row][slice]->data[1]))
                    return false;
            }
        }
        return true;
    }
    
    void set_old(){
        old_list_map = list_map;
        old_key_map = key_map;
        old_size = size;
        set_old_data();
        
        old_key_map = key_map;
        
        sep->set_update(false);
        perm_sep->set_update(false);
    }
    
    
    void revert_old(){
        revert_old_data();
        size = old_size;
        list_map = old_list_map;
        key_map = old_key_map;
        size.resize(list_map.size());
        for(typename std::vector<std::vector<node_ptr> >::iterator it = list_map.begin(); it != list_map.end(); it++)
            size[it-list_map.begin()] = (*it).size();
        sep->reject_update();
        perm_sep->reject_update();
    }
    
    
    void remove_old_list_references(iVector rows_to_remove){
        for(iVector::iterator row = rows_to_remove.begin(); row != rows_to_remove.end(); row++){
            for(int slice = 0; slice < old_size[0]; slice++){
                old_list_map[*row][slice]->right_node.reset();
                old_list_map[*row][slice]->left_node.reset();
                old_list_map[*row][slice]->old_right_node.reset();
                old_list_map[*row][slice]->old_left_node.reset();
            }
        }
    }
    
    std::pair<int, int> get_bead_indices(size_t key){
        auto it = key_map.find(key);
        int row = it->second->row_number;
        int col = it->second->column_number;
        return std::pair<int, int>(row,col);
    }
    
    size_t get_prev_bead_key(size_t key, int dist_back){
        auto it = key_map.find(key);
        int row = it->second->row_number;
        int col = it->second->column_number;
        
        if(dist_back > col){
            col = col-dist_back+size[0];
            row = list_map[row][0]->left_node->row_number;
        }
        else
            col = col - dist_back;
        
        return list_map[row][col]->key;
    }
    
    size_t get_next_bead_key(size_t key, int dist_forward){
        auto it = key_map.find(key);
        int row = it->second->row_number;
        int col = it->second->column_number;
        
        if(dist_forward + col >= size[0]){
            col = (col + dist_forward)%size[0];
            row = list_map[row][size[0]-1]->right_node->row_number;
        }
        else
            col = col+dist_forward;
        return list_map[row][col]->key;
    }
    
    size_t get_bead_key(int row, int col){
        return list_map[row][col]->key;
    }
};

#endif
