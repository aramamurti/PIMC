//
//  beadContainer.h
//  BEC-monopoles
//
//  Created by Adith Ramamurti on 6/8/15.
//  Copyright (c) 2015 Adith Ramamurti. All rights reserved.
//

#ifndef BEC_monopoles_beadContainer_h
#define BEC_monopoles_beadContainer_h
#include "uni_header.h"

template<class T>

class PathList {
private:
    
    class Node {
    public:
        typedef boost::shared_ptr<Node> node_ptr;
        
        int row_number;
        int column_number;
        T data;
        T olddat;
        node_ptr leftNode;
        node_ptr rightNode;
        node_ptr oldRightNode;
        Node(T e) {
            this->data = e;
            leftNode = node_ptr();
            rightNode = node_ptr();
        }
        ~Node() {
        }
    };
    typedef boost::shared_ptr<Node> node_ptr;
    
    bool circular = false;
    std::vector<int> size;
    std::vector<node_ptr> head;
    std::vector<node_ptr> tail;
    std::vector<int> permute_start_orders;
    std::vector<int> permute_end_order;
    
    std::vector<int> pso;
    std::vector<int> peo;
    std::vector<int> rso;
    
    std::vector<int> repermute_start_order;
    std::vector<int> repermute_end_order;
    
    int permute_position;
    
    std::vector<std::vector<node_ptr> > list_map;
    std::vector<std::vector<node_ptr> > old_list_map;
    
    
public:
    
    PathList() {
        size.resize(0);
        head.resize(0);
        tail.resize(0);
    }
    
    PathList(std::string infile)
    {
        size.resize(0);
        head.resize(0);
        tail.resize(0);
        
        std::ifstream in(infile.c_str());
        std::string s;
        size_t numbreaks = 0;
        size_t numline = 0;
        std::vector<std::vector<std::vector<float> > > dataset;
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
                        std::vector<float> data(n);
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
                        std::vector<std::vector<float> > dset2 = dataset[numline];
                        std::vector<float> nodedat= dset2[std::atoi(ent[1].c_str())];
                        pushBack(nodedat, numline, std::atoi(ent[0].c_str()),std::atoi(ent[1].c_str()));
                    }
                }
                token = s.substr(1,s.length()-2);
                if(numbreaks == 0){
                    size_t n = std::count(token.begin(),token.end(),',')+1;
                    std::vector<float> data(n);
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
                    std::vector<std::vector<float> > dset2 = dataset[numline];
                    std::vector<float> nodedat= dset2[std::atoi(ent[1].c_str())];
                    pushBack(nodedat, numline, std::atoi(ent[0].c_str()),std::atoi(ent[1].c_str()));
                }
                
                if(numbreaks == 2){
                    std::vector<std::string> ent;
                    boost::split(ent, s, boost::is_any_of("("), boost::token_compress_on);
                    std::string ss = ent.back();
                    ss = ss.substr(0,s.length()-2);
                    boost::split(ent, ss, boost::is_any_of(", "), boost::token_compress_on);
                    tail[numline]->rightNode = head[std::atoi(ent[0].c_str())];
                    tail[numline]->rightNode->leftNode = tail[numline];
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
        for(typename std::vector<node_ptr>::iterator it = head.begin(); it != head.end(); it++){
            node_ptr temp = *it;
            while (temp != NULL) {
                node_ptr next(temp->rightNode);
                temp->rightNode.reset();
                temp->leftNode.reset();
                temp.reset();
                temp = next;
            }
        }
    }
    
    void pushBack(T t, int index = 0, int row_number = -1, int column_number = -1) {
        node_ptr Node(new class Node(t));
        if(row_number == -1)
            Node->row_number = index;
        else
            Node->row_number = row_number;
        
        if(size.size() <= index){
            size.resize(index+1);
            head.resize(index+1);
            tail.resize(index+1);
        }
        ++size[index];
        
        list_map.resize(size.size());
        list_map[index].resize(size[index]);
        
        if(column_number == -1)
            Node->column_number = size[index]-1;
        else
            Node->column_number = column_number;
        
        if (size[index] == 1) {
            head[index] = tail[index] = Node;
            list_map[index][size[index]-1] = Node;
            return;
        }
        Node->leftNode = tail[index];
        tail[index]->rightNode = Node;
        Node->leftNode = tail[index];
        tail[index] = Node;
        
        list_map[index][size[index]-1] = Node;
    }
    
    void shift_all(int row, T shift){
        node_ptr temp = head[row];
        for(typename T::iterator it = temp->data.begin(); it != temp->data.end(); it++){
            *it += shift[it-temp->data.begin()];
        }
        temp = temp->rightNode;
        while(temp != head[row]){
            for(typename T::iterator it = temp->data.begin(); it != temp->data.end(); it++){
                *it += shift[it-temp->data.begin()];
            }
            temp = temp->rightNode;
        }
    }
    
    void set_old_data(){
        for(int i = 0; i < size.size(); i++){
            for(int j = 0; j<size[i];j++){
                list_map[i][j]->olddat = list_map[i][j]->data;
            }
        }
    }
    
    void revert_old_data(){
        for(int i = 0; i < size.size(); i++){
            for(int j = 0; j<size[i];j++){
                list_map[i][j]->data = list_map[i][j]->olddat;
            }
        }
    }
    
    void set_permutation(std::vector<int> i, std::vector<int> j, int pos, int dist = 0){
        permute_position = (pos+dist);
        permute_start_orders.resize(i.size());
        permute_end_order.resize(j.size());
        
        std::copy(i.begin(), i.end(), permute_start_orders.begin());
        std::copy(j.begin(), j.end(), permute_end_order.begin());
        
        if(permute_position > size[0]){
            std::vector<std::vector<int> > perms = circular_perm(i,j);
            repermute_start_order = perms[0];
            repermute_end_order = perms[1];
        }
        else{
            repermute_start_order = i;
            repermute_end_order = j;
        }
    }
    
    void permute(bool reverse = false){
        std::vector<int> i;
        std::vector<int> j;
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
                    *it = temps[j[it-temps2.begin()]]->rightNode;
                
                for(typename std::vector<int>::iterator it = i.begin(); it != i.end(); it++){
                    if(permute_position != 0){
                        temps[*it]->oldRightNode = temps[*it]->rightNode;
                        temps[*it]->rightNode = temps2[it-i.begin()];
                        temps[*it]->rightNode->leftNode = temps[*it];
                    }
                }
                
                for(int curpos = pos; curpos < size[0]; curpos++){
                    for(std::vector<int>::iterator it = i.begin(); it != i.end();it++){
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
                
                for(typename std::vector<int>::iterator it = i.begin(); it != i.end(); it++){
                    if(permute_position != 0){
                        temps[*it]->rightNode = temps[*it]->oldRightNode;
                        temps[*it]->rightNode->leftNode = temps[*it];
                    }
                }
                list_map = old_list_map;
            }
        }
    }
    
    std::vector<std::vector<int> > circular_perm(std::vector<int> lc, std::vector<int> end){
        
        std::vector<node_ptr> temps;
        
        for(int ptc = 0; ptc < size.size(); ptc ++)
            temps.push_back(list_map[ptc][size[ptc]-1]);
        
        std::vector<int> nlc;
        std::vector<int> nend;
        for(std::vector<int>::iterator it = lc.begin(); it != lc.end(); it++){
            nlc.push_back(temps[*it]->rightNode->row_number);
            nend.push_back(temps[end[it-lc.begin()]]->rightNode->row_number);
        }
        
        std::vector<std::vector<int> > reperms;
        reperms.push_back(nlc);
        reperms.push_back(nend);
        return reperms;
    }
    
    std::vector<float> get_bead_data(int row, int slice){
        
        if(slice < size[row]){
            return list_map[row][slice]->data;
        }
        else if(circular){
            slice = slice%size[row];
            row = list_map[row][size[row]-1]->rightNode->row_number;
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
            row = list_map[row][size[row]-1]->rightNode->row_number;
            list_map[row][slice]->data = data;
        }
    }
    
    std::vector<T> get_pair_same_path(int row, int start, int dist){
        
        int row2 = row;
        int end = start+dist;
        
        if(start >= size[row] && end >= size[row] && circular){
            start = start%size[row];
            row = list_map[row][size[row]-1]->rightNode->row_number;
        }

        if(end >= size[row] && circular){
            end = end%size[row];
            row2 = list_map[row2][size[row]-1]->rightNode->row_number;
        }
        
        std::vector<T> ret(0);
        ret.push_back(list_map[row][start]->data);
        ret.push_back(list_map[row2][end]->data);
        return ret;
    }
    
    std::vector<T> get_pair_same_slice(int row1, int row2, int slice){
        
        if(slice >= size[row1] && circular){
            slice = slice%size[row1];
            row1 = list_map[row1][size[row1]-1]->rightNode->row_number;
            row2 = list_map[row2][size[row2]-1]->rightNode->row_number;
        }
        std::vector<T> ret(0);
        ret.push_back(list_map[row1][slice]->data);
        ret.push_back(list_map[row2][slice]->data);
        return ret;
        
    }
    
    std::vector<int> get_cycles(){
        std::vector<int> cyclenum(size.size(),0);
        for(int row = 0; row < size.size(); row++){
            node_ptr temp = head[row]->rightNode;
            int tempPos = 1;
            while(temp != head[row]){
                temp = temp->rightNode;
                ++tempPos;
            }
            int cycNum = tempPos/size[row];
            
            cyclenum[cycNum-1]++;
        }
        
        return cyclenum;
    }
    
    
    
    void make_circular(bool circ = true){
        
        if(circ){
            for(typename std::vector<node_ptr>::iterator it = head.begin(); it != head.end(); it++){
                node_ptr ln = *it;
                ln->leftNode = tail[it-head.begin()];
                tail[it-head.begin()]->rightNode = ln;
            }
            circular = true;
        }
        else{
            for(typename std::vector<node_ptr>::iterator it = head.begin(); it != head.end(); it++){
                node_ptr ln = *it;
                ln->leftNode.reset();
                tail[it-head.begin()]->rightNode.reset();
            }
            circular = false;
        }
    }
    
    
    
    void print_list(int index = 0){
        node_ptr temp = head[index];
        if(circular){
            bool headprint = false;
            int numPrint = 0;
            while (!headprint || numPrint < size[index]) {
                if(temp == head[index])
                    headprint = true;
                T beep = temp->data;
                for(typename T::iterator it = beep.begin(); it != beep.end(); it++){
                    if(it == beep.begin())
                        std::cout<<"(";
                    std::cout << *it;
                    if(it != beep.end()-1)
                        std::cout <<", ";
                    
                }
                std::cout <<")";
                temp = temp->rightNode;
                numPrint++;
                if(numPrint < size[index])
                    std::cout<< ", ";
            }
            
        }
        else{
            while (temp != NULL) {
                T beep = temp->data;
                for(typename T::iterator it = beep.begin(); it != beep.end(); it++)
                    std::cout << *it << ", ";
                temp = temp->rightNode;
            }
        }
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
                std::cout << ", ";
            
        }
        std::cout<<std::endl;
    }
    
    void print_list_indices(int index = 0){
        node_ptr temp = head[index];
        if(circular){
            bool headprint = false;
            int numPrint = 0;
            while (!headprint || numPrint < size[index]) {
                if(temp == head[index])
                    headprint = true;
                std::cout << "("<<temp->row_number<<", "<<temp->column_number<<")";
                temp = temp->rightNode;
                numPrint++;
                if(numPrint < size[index])
                    std::cout<< ", ";
            }
            
        }
        else{
            while (temp != NULL) {
                std::cout << "("<<temp->row_number<<", "<<temp->column_number<<"), " << ", ";
                temp = temp->rightNode;
            }
        }
        std::cout << std::endl;
    }
    
    void print_list_next_indices(int index = 0){
        node_ptr temp = head[index];
        if(circular){
            bool headprint = false;
            int numPrint = 0;
            while (!headprint || numPrint < size[index]) {
                if(temp == head[index])
                    headprint = true;
                std::cout << "("<<temp->rightNode->row_number<<", "<<temp->rightNode->column_number<<")";
                temp = temp->rightNode;
                numPrint++;
                if(numPrint < size[index])
                    std::cout<< ", ";
            }
            
        }
        else{
            while (temp != NULL) {
                std::cout << "("<<temp->rightNode->row_number<<", "<<temp->rightNode->column_number<<"), " << ", ";
                temp = temp->rightNode;
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
    
    
    std::vector<int> get_rep_swap_order(){
        return rso;
    }
    
    void set_prev_perm(){
        pso = permute_start_orders;
        peo = permute_end_order;
        rso = repermute_start_order;
    }
    
    void reset_permute(){
        permute_start_orders = pso;
        permute_end_order = peo;
        repermute_start_order = rso;
    }
};

    
#endif
