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
        size_t key;
        
        T data;
        T olddat;
        
        node_ptr left_node;
        node_ptr right_node;
        node_ptr old_right_node;
        
        Node(T e) {
            this->data = e;
            left_node = node_ptr();
            right_node = node_ptr();
        }
        ~Node() {
        }
    };
    
    typedef boost::shared_ptr<Node> node_ptr;
    boost::hash<boost::tuple<int, int, int> > triple_hash;
    
    bool circular = false;
    
    int push_counter = 0;
    iVector size;
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
    
    std::vector<std::vector<node_ptr> > worm;
    std::vector<std::vector<node_ptr> > ye_olde_worm;
    
    node_ptr worm_head;
    node_ptr worm_tail;
    int worm_size = 0;
    std::pair<int, int> worm_head_index;
    std::pair<int, int> worm_tail_index;
    
    node_ptr old_worm_head;
    node_ptr old_worm_tail;
    int old_worm_size = 0;
    std::pair<int, int> old_worm_head_index;
    std::pair<int, int> old_worm_tail_index;
    
    boost::unordered_map<size_t, node_ptr> key_map;
    
    boost::shared_ptr<Separation_Table> sep;
    boost::shared_ptr<Neighbor_Table<node_ptr> > nt;
    
public:
    
    PathList(int ndim, double boxsize = -1) {
        size.resize(0);
        sep = boost::shared_ptr<Separation_Table>(new Separation_Table(boxsize));
        nt = boost::shared_ptr<Neighbor_Table<node_ptr> >(new Neighbor_Table<node_ptr>(boxsize, ndim));
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
        else
            node->row_number = row_number;
        
        if(size.size() <= index){
            size.resize(index+1);
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
            nt->add_bead(node);
        }
        else
            node->key = key;
        
        update_key_map(node);
        
        
        if (size[index] == 1) {
            list_map[index][size[index]-1] = node;
            return;
        }
        
        node->left_node = list_map[index][size[index]-2];
        node->left_node->right_node = node;
        
        list_map[index][size[index]-1] = node;
    }
    
    void set_key_add_sep(node_ptr node){
        boost::tuple<int, int,int> rc(push_counter, node->row_number, node->column_number);
        node->key = triple_hash(rc);
        
        push_counter++;
        sep->add_bead(node->column_number, node->key, node->data);
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
        if(it != key_map.end())
            key_map.erase(it);
        key_map.insert(std::pair<size_t, node_ptr>(key, node));
        
    }
    
    void generate_separations(){
        sep->initialize_sep_table();
    }
    
    void generate_neighbors(){
        nt->set_up_neighbor_table();
    }
    
    void shift_all(int row, T shift){
        node_ptr temp = list_map[row][0];
        std::vector<std::pair<size_t, int> > key_list;
        std::vector<T> new_locs;
        
        for(typename T::iterator it = temp->data.begin(); it != temp->data.end(); it++){
            *it += shift[it-temp->data.begin()];
        }
        nt->update_bead(temp);
        
        new_locs.push_back(temp->data);
        key_list.push_back(std::pair<size_t, int>(temp->key, temp->column_number));
        
        temp = temp->right_node;
        while(temp != list_map[row][0]){
            for(typename T::iterator it = temp->data.begin(); it != temp->data.end(); it++){
                *it += shift[it-temp->data.begin()];
            }
            
            nt->update_bead(temp);
            
            new_locs.push_back(temp->data);
            key_list.push_back(std::pair<size_t, int>(temp->key, temp->column_number));
            
            temp = temp->right_node;
        }
        
        for(std::vector<std::pair<size_t, int> >::iterator it = key_list.begin(); it!= key_list.end(); it++){
            sep->update_bead(it->second, it->first, new_locs[it-key_list.begin()]);
        }
    }
    
    void set_old_data(){
        for(int i = 0; i < size.size(); i++){
            for(int j = 0; j<size[i];j++){
                list_map[i][j]->olddat = list_map[i][j]->data;
            }
        }
        sep->set_old_maps();
        nt->set_old_table();
    }
    
    void revert_old_data(){
        for(int i = 0; i < size.size(); i++){
            for(int j = 0; j<size[i];j++){
                list_map[i][j]->data = list_map[i][j]->olddat;
            }
        }
        sep->revert_maps();
        nt->reset_old_table();
    }
    
    void set_permutation(iVector i, iVector j, int pos, int dist = 0){
        permute_position = (pos+dist);
        permute_start_orders.resize(i.size());
        permute_end_order.resize(j.size());
        
        std::copy(i.begin(), i.end(), permute_start_orders.begin());
        std::copy(j.begin(), j.end(), permute_end_order.begin());
        
        if(permute_position > size[0]){
            iiVector perms = circular_perm(i,j);
            repermute_start_order = perms[0];
            repermute_end_order = perms[1];
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
            }
        }
    }
    
    iiVector circular_perm(iVector lc, iVector end){
        
        std::vector<node_ptr> temps;
        
        for(int ptc = 0; ptc < size.size(); ptc ++)
            temps.push_back(list_map[ptc][size[ptc]-1]);
        
        iVector nlc;
        iVector nend;
        for(iVector::iterator it = lc.begin(); it != lc.end(); it++){
            nlc.push_back(temps[*it]->right_node->row_number);
            nend.push_back(temps[end[it-lc.begin()]]->right_node->row_number);
        }
        
        iiVector reperms;
        reperms.push_back(nlc);
        reperms.push_back(nend);
        return reperms;
    }
    
    dVector get_bead_data(int row, int slice){
        
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
    
    void set_bead_data(int row, int slice, T data){
        if(slice < size[row]){
            list_map[row][slice]->data = data;
        }
        else if(circular){
            slice = slice%size[row];
            row = list_map[row][size[row]-1]->right_node->row_number;
            list_map[row][slice]->data = data;
        }
        sep->update_bead(slice, list_map[row][slice]->key, data);
        nt->update_bead(list_map[row][slice]);
    }
    
    std::vector<T> get_pair_same_path(int row, int start, int dist){
        
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
        
        std::vector<T> ret(0);
        ret.push_back(list_map[row][start]->data);
        ret.push_back(list_map[row2][end]->data);
        return ret;
    }
    
    std::vector<T> get_pair_same_slice(int row1, int row2, int slice){
        
        if(slice >= size[row1] && circular){
            slice = slice%size[row1];
            row1 = list_map[row1][size[row1]-1]->right_node->row_number;
            row2 = list_map[row2][size[row2]-1]->right_node->row_number;
        }
        std::vector<T> ret(0);
        ret.push_back(list_map[row1][slice]->data);
        ret.push_back(list_map[row2][slice]->data);
        return ret;
        
    }
    
    T get_path_separation(int row1, int row2, int slice){
        return sep->get_separation(std::pair<size_t,size_t>(list_map[row1][slice]->key,list_map[row2][slice]->key));
    }
    
    T get_worm_separation(int worm_row1, int worm_row2, int slice){
        if(worm_row1 == worm_head_index.first || worm_row2 == worm_head_index.first){
            if(slice < worm_head_index.second)
                return T(0);
        }
        if(worm_row1 == worm_tail_index.first || worm_row2 == worm_tail_index.first){
            if(slice > worm_tail_index.second)
                return T(0);
        }
        return sep->get_separation(std::pair<size_t,size_t>(worm[worm_row1][slice]->key,worm[worm_row2][slice]->key));
    }
    
    T get_worm_path_separation(int path_row, int worm_row, int slice){
        if(worm_row == worm_head_index.first){
            if(slice < worm_head_index.second)
                return T(0);
        }
        if(worm_row == worm_tail_index.first){
            if(slice > worm_tail_index.second)
                return T(0);
        }
        
        return sep->get_separation(std::pair<size_t,size_t>(list_map[path_row][slice]->key,worm[worm_row][slice]->key));
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
    
    
    /******************
     Worm methods
     ******************/
    
    void initialize_worm(){
        worm.resize(1);
        worm[0].resize(size[0]);
        worm_size = 0;
    }
    
    void new_worm(T t, int slice, size_t key = -1){
        if(worm_size == 0){
            node_ptr node(new Node(t));
            node->row_number = 0;
            node->column_number = slice;
            worm_head_index.second = worm_tail_index.second = slice;
            worm_head_index.first = worm_tail_index.first = 0;
            worm_head = worm_tail = node;
            worm[0][slice] = worm_head;
            if(key == -1){
                set_key_add_sep(node);
                nt->add_bead(node);
            }
            else
                node->key = key;
            
            update_key_map(node);
            
            worm_size++;
        }
    }
    
    void worm_push_back(T t, size_t key = -1){
        node_ptr node(new Node(t));
        if(worm_size == 0){
            std::cout << "No worm exists! Can't push back!" << std::endl;
            return;
        }
        else{
            worm_tail->right_node = node;
            node->left_node = worm_tail;
            worm_tail = node;
            
            worm_tail_index.second++;
            if(worm_tail_index.second >= worm[0].size()){
                worm_tail_index.first++;
                worm_tail_index.second = 0;
            }
            if(worm_tail_index.first >= worm.size()){
                worm.resize(worm_tail_index.first+1);
                worm[worm_tail_index.first].resize(worm[0].size());
            }
            worm[worm_tail_index.first][worm_tail_index.second] = worm_tail;
            worm_tail->row_number = worm_tail_index.first;
            worm_tail->column_number = worm_tail_index.second;
            if(key == -1){
                set_key_add_sep(node);
                nt->add_bead(node);
            }
            else
                node->key = key;
            
            update_key_map(node);
            
            
            worm_size++;
        }
    }
    
    void worm_push_front(T t, size_t key = -1){
        node_ptr node(new Node(t));
        if(worm_size == 0){
            std::cout << "No worm exists! Can't push front!" << std::endl;
            return;
        }
        else{
            worm_head->left_node = node;
            node->right_node = worm_head;
            worm_head = node;
            
            worm_head_index.second--;
            if(worm_head_index.second < 0){
                worm_tail_index.first++;
                worm_head_index.second = worm[0].size()-1;
                worm.insert(worm.begin(),std::vector<node_ptr>(worm[0].size()));
                worm[worm_head_index.first][worm_head_index.second] = node;
                node_ptr temp_node = node->right_node;
                for(int i = 0; i < worm_size; i++){
                    temp_node->row_number++;
                    temp_node = temp_node->right_node;
                }
            }
            else{
                worm[worm_head_index.first][worm_head_index.second] = node;
            }
            
            worm_head->row_number = worm_head_index.first;
            worm_head->column_number = worm_head_index.second;
            if(key == -1){
                set_key_add_sep(node);
                nt->add_bead(node);
            }
            else
                node->key = key;
            
            update_key_map(node);
            
            worm_size++;
        }
    }
    
    T worm_pop_back(bool remove = false){
        if(worm_size == 0){
            std::cout << "No worm exists! Can't pop back!" << std::endl;
            return T(0);
        }
        else if(worm_size == 1){
            node_ptr node = worm_head;
            worm_head.reset();
            worm_tail.reset();
            worm[worm_tail_index.first][worm_tail_index.second] = NULL;
            worm_tail_index.first = 0;
            worm_head_index.second = 0;
            worm_tail_index.second = 0;
            worm_size = 0;
            worm.resize(0);
            if(remove)
                sep->remove_bead(node->column_number, node->key);
            return node->data;
        }
        else{
            node_ptr node = worm_tail;
            worm_tail = worm_tail->left_node;
            worm_tail->right_node.reset();
            node->left_node.reset();
            worm[worm_tail_index.first][worm_tail_index.second] = NULL;
            worm_tail_index.second--;
            if(worm_tail_index.second < 0){
                worm_tail_index.second = worm[worm_tail_index.first-1].size()-1;
                worm_tail_index.first--;
                worm.erase(worm.end());
            }
            worm_size--;
            if(remove){
                sep->remove_bead(node->column_number, node->key);
                nt->remove_bead(node);
            }
            return node->data;
        }
    }
    
    T worm_pop_front(bool remove = false){
        if(worm_size == 0){
            std::cout << "No worm exists! Can't pop front!" << std::endl;
            return T(0);
        }
        else if(worm_size == 1){
            node_ptr node = worm_head;
            worm_head.reset();
            worm_tail.reset();
            worm[worm_tail_index.first][worm_tail_index.second] = NULL;
            worm_tail_index.first = 0;
            worm_tail_index.second = 0;
            worm_head_index.second = 0;
            worm_size = 0;
            worm.resize(0);
            if(remove)
                sep->remove_bead(node->column_number, node->key);
            return node->data;
        }
        else{
            node_ptr node = worm_head;
            worm_head = worm_head->right_node;
            worm_head->left_node.reset();
            node->right_node.reset();
            worm[worm_head_index.first][worm_head_index.second] = NULL;
            worm_head_index.second++;
            if(worm_head_index.second >= worm[0].size()){
                worm_head_index.second = 0;
                worm_tail_index.first--;
                worm.erase(worm.begin());
                node_ptr temp_node = worm_head;
                for(int i = 0; i < worm_size-1; i++){
                    temp_node->row_number--;
                    temp_node = temp_node->right_node;
                }
            }
            worm_size--;
            if(remove){
                sep->remove_bead(node->column_number, node->key);
                nt->remove_bead(node);
            }
            return node->data;
        }
    }
    
    void move_path_to_worm(int row, int head_col, int tail_dist_back){
        set_old_worm();
        set_old_list();
        
        if(worm_size != 0){
            std::cout << "Worm exists! Cannot move a path into the worm." << std::endl;
            return;
        }
        else{
            initialize_worm();
            int col_pos = head_col-1;
            new_worm(list_map[row][head_col]->data, head_col,list_map[row][head_col]->key);
            node_ptr temp_node = list_map[row][head_col]->right_node;
            std::vector<int> rows_to_remove(1,row);
            while(temp_node != list_map[row][head_col]){
                if(col_pos == size[row]){
                    col_pos = 0;
                    rows_to_remove.push_back(temp_node->row_number);
                }
                worm_push_back(temp_node->data, temp_node->key);
                temp_node = temp_node->right_node;
                col_pos++;
            }
            
            while(tail_dist_back != 0){
                worm_pop_back(true);
                tail_dist_back--;
            }
            std::sort(rows_to_remove.begin(), rows_to_remove.end());
            auto last = std::unique(rows_to_remove.begin(), rows_to_remove.end());
            rows_to_remove.erase(last, rows_to_remove.end());
            
            for(std::vector<int>::iterator it = rows_to_remove.begin(); it != rows_to_remove.end(); it++){
                *it = *it - (it-rows_to_remove.begin());
                list_map.erase(list_map.begin()+ *it);
                size.erase(size.begin()+* it);
            }
            
            reset_indices(rows_to_remove);
        }
    }
    
    void move_worm_to_path(){
        set_old_worm();
        set_old_list();
        
        if(worm_size == 0){
            std::cout << "No worm exists." << std::endl;
            return;
        }
        else{
            int first_row = size.size();
            int row = size.size();
            
            bool head_tail_same_row;
            if(worm_tail_index.second>=worm_head_index.second){
                while(worm_head_index.second != 0)
                    worm_push_front(T(0));
                while(worm_tail_index.second != worm[0].size()-1)
                    worm_push_back(T(0));
                int colpos = 0;
                int maxcol = size[0]-1;
                node_ptr temp = worm[0][0];
                while(temp != NULL){
                    if(colpos > maxcol){
                        row++;
                        colpos = 0;
                    }
                    push_back(temp->data, row, temp->key);
                    temp = temp->right_node;
                    colpos++;
                    worm_pop_front();
                }
            }
            else{
                int difference = worm_head_index.second - worm_tail_index.second;
                difference--;
                while(difference != 0){
                    worm_push_front(T(0));
                    difference--;
                }
                int maxcol = size[0]-1;
                node_ptr temp = worm_head;
                int colpos = worm_head_index.second-1;
                for(int i = 0; i <= colpos; i++){
                    push_back(worm[worm_tail_index.first][i]->data, row, worm[worm_tail_index.first][i]->key);
                }
                colpos++;
                while(temp != NULL & temp->row_number != worm_tail_index.first){
                    if(colpos > maxcol){
                        row++;
                        colpos = 0;
                    }
                    push_back(temp->data, row, temp->key);
                    temp = temp->right_node;
                    colpos++;
                    worm_pop_front();
                }
                while(temp != NULL){
                    temp = temp->right_node;
                    worm_pop_front();
                }
            }
            
            list_map[first_row][0]->left_node = list_map[row][size[row]-1];
            list_map[first_row][0]->left_node->right_node = list_map[first_row][0];
            for(row = first_row; row < size.size()-1; row++){
                list_map[row][size[row]-1]->right_node = list_map[row+1][0];
                list_map[row][size[row]-1]->right_node->left_node = list_map[row][size[row]-1];
            }
        }
    }
    
    void reset_indices(iVector rtr){
        std::vector<int> new_indices(size.size()+rtr.size(),0);
        for(std::vector<int>::iterator it2 = new_indices.begin(); it2 != new_indices.end(); it2 ++)
            *it2 = it2 - new_indices.begin();
        
        for(std::vector<int>::iterator it1= rtr.begin(); it1 != rtr.end(); it1++)
            for(std::vector<int>::iterator it2 = new_indices.begin(); it2 != new_indices.end(); it2 ++)
                if(*it2  > *it1)
                    *it2 = *it2 - 1;
        
        for(typename std::vector<std::vector<node_ptr> >::iterator it = list_map.begin(); it != list_map.end(); it++)
            for(typename std::vector<node_ptr>::iterator it2 = (*it).begin(); it2 != (*it).end(); it2++)
                (*it2)->row_number = new_indices[(*it2)->row_number];
    }
    
    void set_old_worm(){
        ye_olde_worm = worm;
        old_worm_size = worm_size;
        old_worm_head = worm_head;
        old_worm_tail = worm_tail;
        old_worm_head_index = worm_head_index;
        old_worm_tail_index = worm_tail_index;
    }
    
    void set_old_list(){
        old_list_map = list_map;
    }
    
    void set_old_worm_data(){
        node_ptr temp = worm[worm_head_index.first][worm_head_index.second];
        while(temp != NULL){
            temp->olddat = temp->data;
            temp = temp->right_node;
        }
        sep->set_old_maps();
        nt->set_old_table();
    }
    
    void revert_worm_old_data(){
        node_ptr temp = worm[worm_head_index.first][worm_head_index.second];
        while(temp != NULL){
            temp->data = temp->olddat;
            temp = temp->right_node;
        }
        sep->revert_maps();
        nt->reset_old_table();
    }
    
    void reset_worm(){
        worm = ye_olde_worm;
        worm_size = old_worm_size;
        worm_head = old_worm_head;
        worm_tail = old_worm_tail;
        worm_head_index = old_worm_head_index;
        worm_tail_index = old_worm_tail_index;
    }
    
    void reset_list(){
        list_map = old_list_map;
        size.resize(list_map.size());
        for(typename std::vector<std::vector<node_ptr> >::iterator it = list_map.begin(); it != list_map.end(); it++)
            size[it-list_map.begin()] = (*it).size();
    }
    
    void swap_into_worm_head(int row, int dist){
        if(worm_size == 0)
            return;
        else{
            int swap_bead_index = (worm_head_index.second-dist+size[row])%size[row];
            for(int i = 0; i < dist - 1; i++)
                worm_push_front(T(0));
            worm_push_front(list_map[row][swap_bead_index]->data, list_map[row][swap_bead_index]->key);
            int col_pos = (worm_tail_index.second-dist+size[row])%size[row];
            node_ptr temp_node = list_map[row][swap_bead_index]->left_node;
            std::vector<int> rows_to_remove(1,row);
            while(temp_node != list_map[row][swap_bead_index]){
                if(col_pos == 0){
                    col_pos = size[row];
                    rows_to_remove.push_back(temp_node->row_number);
                }
                worm_push_front(temp_node->data, temp_node->key);
                temp_node = temp_node->left_node;
                col_pos--;
            }
            
            dist--;
            while(dist != 0){
                worm_pop_front(true);
                dist--;
            }
            std::sort(rows_to_remove.begin(), rows_to_remove.end());
            auto last = std::unique(rows_to_remove.begin(), rows_to_remove.end());
            rows_to_remove.erase(last, rows_to_remove.end());
            
            for(std::vector<int>::iterator it = rows_to_remove.begin(); it != rows_to_remove.end(); it++){
                *it = *it - (it-rows_to_remove.begin());
                list_map.erase(list_map.begin()+ *it);
                size.erase(size.begin()+* it);
            }
            reset_indices(rows_to_remove);
            
        }
    }
    
    void swap_into_worm_tail(int row, int dist){
        if(worm_size == 0)
            return;
        else{
            int swap_bead_index = (worm_tail_index.second+dist)%size[row];
            for(int i = 0; i < dist - 1; i++)
                worm_push_back(T(0));
            worm_push_back(list_map[row][swap_bead_index]->data, list_map[row][swap_bead_index]->key);
            int col_pos = (worm_tail_index.second+dist)%size[row];
            node_ptr temp_node = list_map[row][swap_bead_index]->right_node;
            std::vector<int> rows_to_remove(1,row);
            while(temp_node != list_map[row][swap_bead_index]){
                if(col_pos == size[row]){
                    col_pos = 0;
                    rows_to_remove.push_back(temp_node->row_number);
                }
                worm_push_back(temp_node->data, temp_node->key);
                temp_node = temp_node->right_node;
                col_pos++;
            }
            
            dist--;
            while(dist != 0){
                worm_pop_back(true);
                dist--;
            }
            std::sort(rows_to_remove.begin(), rows_to_remove.end());
            auto last = std::unique(rows_to_remove.begin(), rows_to_remove.end());
            rows_to_remove.erase(last, rows_to_remove.end());
            
            for(std::vector<int>::iterator it = rows_to_remove.begin(); it != rows_to_remove.end(); it++){
                *it = *it - (it-rows_to_remove.begin());
                list_map.erase(list_map.begin()+ *it);
                size.erase(size.begin()+ *it);
            }
            
            reset_indices(rows_to_remove);
            
        }
    }
    
    dVector get_worm_bead_data(int row, int slice){
        
        if(slice < worm[row].size() && slice >= worm_head_index.second){
            return worm[row][slice]->data;
        }
        else if(row < worm.size()-2){
            slice = slice%worm[row].size();
            row += 1;
            return worm[row][slice]->data;
        }
        else if(row == worm.size()-2 && slice%worm[row].size() <= worm_tail_index.second){
            slice = slice%worm[row].size();
            row += 1;
            return worm[row][slice]->data;
        }
        else
            return T(0);
    }
    
    void set_worm_bead_data(int row, int slice, T data){
        
        if(slice < worm[row].size() && slice >= worm_head_index.second){
            worm[row][slice]->data = data;
        }
        else if(row < worm.size()-2){
            slice = slice%worm[row].size();
            row += 1;
            worm[row][slice]->data = data;
        }
        else if(row == worm.size()-2 && slice%worm[row].size() <= worm_tail_index.second){
            slice = slice%worm[row].size();
            row += 1;
            worm[row][slice]->data = data;
        }
        else
            return;
        sep->update_bead(slice, worm[row][slice]->key, data);
    }
    
    
    /************
     Printing Methods
     ***********/
    
    void print_separation(int row1, int row2, int slice){
        dVector dist = sep->get_separation(std::pair<size_t,size_t>(list_map[row1][slice]->key,list_map[row2][slice]->key));
        for(dVector::iterator it = dist.begin(); it != dist.end(); it++)
            std::cout << *it <<"\t";
        std::cout<<std::endl;
    }
    
    void print_list_map(int index = 0){
        std::vector<node_ptr> list = list_map[index];
        for(typename std::vector<node_ptr>::iterator it = list.begin(); it != list.end(); it++){
            node_ptr node = *it;
            T data = node->data;
            for(typename T::iterator it2 = data.begin(); it2 != data.end(); it2++){
                if(it2 == data.begin())
                    std::cout<<"(";
                std::cout << *it2;
                if(it2 != data.end()-1)
                    std::cout <<", ";
                
            }
            std::cout <<")";
            if(it != list.end()-1)
                std::cout << ",\t";
            
        }
        std::cout<<std::endl;
    }
    
    void print_list_keys(int index = 0){
        std::vector<node_ptr> list = list_map[index];
        for(typename std::vector<node_ptr>::iterator it = list.begin(); it != list.end(); it++){
            node_ptr node = *it;
            std::cout << node->key;
            if(it != list.end()-1)
                std::cout << ",\t";
        }
        std::cout<<std::endl;
    }
    
    void print_list_indices(int index = 0){
        node_ptr temp = list_map[index][0];
        if(circular){
            bool headprint = false;
            int numPrint = 0;
            while (!headprint || numPrint < size[index]) {
                if(temp == list_map[index][0])
                    headprint = true;
                std::cout << "("<<temp->row_number<<", "<<temp->column_number<<")";
                temp = temp->right_node;
                numPrint++;
                if(numPrint < size[index])
                    std::cout<< ",\t";
            }
            
        }
        else{
            while (temp != NULL) {
                std::cout << "("<<temp->row_number<<", "<<temp->column_number<<"), " << ", ";
                temp = temp->right_node;
            }
        }
        std::cout << std::endl;
        
    }
    void print_worm_data(){
        if(worm_size == 0){
            std::cout<< "Worm is empty!" << std::endl;
            return;
        }
        else
            for(int row = 0; row <= worm_tail_index.first; row++){
                for(int col = 0; col < worm[0].size(); col++){
                    if(worm[row][col] == NULL){
                        std::cout << "NULL";
                    }
                    else{
                        T data = worm[row][col]->data;
                        for(typename T::iterator it2 = data.begin(); it2 != data.end(); it2++){
                            if(it2 == data.begin())
                                std::cout<<"(";
                            std::cout << *it2;
                            if(it2 != data.end()-1)
                                std::cout <<", ";
                            
                        }
                        std::cout <<")";
                    }
                    if(col != worm[0].size() -1)
                        std::cout << ",\t";
                }
                std::cout << std::endl;
            }
    }
    
    void print_list_next_indices(int index = 0){
        node_ptr temp = list_map[index][0];
        if(circular){
            bool headprint = false;
            int numPrint = 0;
            while (!headprint || numPrint < size[index]) {
                if(temp == list_map[index][0])
                    headprint = true;
                std::cout << "("<<temp->right_node->row_number<<", "<<temp->right_node->column_number<<")";
                temp = temp->right_node;
                numPrint++;
                if(numPrint < size[index])
                    std::cout<< ",\t";
            }
            
        }
        else{
            while (temp != NULL) {
                std::cout << "("<<temp->right_node->row_number<<", "<<temp->right_node->column_number<<"), " << ", ";
                temp = temp->right_node;
            }
        }
        std::cout << std::endl;
    }
    
    void print_list_file(int step){
        int numrows = size.size();
        std::stringstream sstm;
        sstm << "config_step_" << step <<".txt";
        std::string result = sstm.str();
        
        std::ofstream out(result.c_str());
        std::streambuf *coutbuf = std::cout.rdbuf();
        std::cout.rdbuf(out.rdbuf());
        for(int j = 0; j < numrows; j++)
            print_list_map(j);
        std::cout << "\n";
        for(int j = 0; j < numrows; j++)
            print_list_indices(j);
        std::cout << "\n";
        for(int j = 0; j < numrows; j++)
            print_list_next_indices(j);
        std::cout.rdbuf(coutbuf);
    }
    
    iVector get_size(){
        return size;
    }
    
    iVector get_worm_dims(){
        iVector worm_dims;
        worm_dims.push_back(worm.size());
        worm_dims.push_back(worm[0].size());
    }
    
    int get_worm_size(){
        return worm_size;
    }
    
    void check_print_neighbor_table(){
        for(int row = 0; row < size.size(); row ++){
            for(int col = 0; col < size[row]; col ++){
                if(nt->check_bead(list_map[row][col])){
                    std::cout << nt->get_bead_grid_num(list_map[row][col]) << ", ";
                }
            }
            std::cout << std::endl;
        }
        
        if(worm_size == 0){
            std::cout<< "Worm is empty!" << std::endl;
            return;
        }
        else
            for(int row = 0; row <= worm_tail_index.first; row++){
                for(int col = 0; col < worm[0].size(); col++){
                    if(worm[row][col] == NULL){
                        std::cout << "NULL";
                    }
                    else{
                        std::cout << nt->get_bead_grid_num(worm[row][col]);
                    }
                    if(col != worm[0].size() -1)
                        std::cout << ", ";
                }
                std::cout << std::endl;
            }
        
    }
};

#endif
