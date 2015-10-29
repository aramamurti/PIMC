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
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>

template<class T>


class LinkedList {
private:
    
    class LinkNode {
    public:
        typedef std::tr1::shared_ptr<LinkNode> node_ptr;
        
        int rownum;
        int colnum;
        T data;
        T olddat;
        node_ptr leftNode;
        node_ptr rightNode;
        node_ptr oldRightNode;
        LinkNode(T e) {
            this->data = e;
            leftNode = node_ptr();
            rightNode = node_ptr();
        }
        ~LinkNode() {
        }
    };
    typedef std::tr1::shared_ptr<LinkNode> node_ptr;
    
    bool circular = false;
    std::vector<int> size;
    std::vector<node_ptr> head;
    std::vector<node_ptr> tail;
    std::vector<int> swapStartOrders;
    std::vector<int> swapEndOrder;
    
    std::vector<int> pso;
    std::vector<int> peo;
    std::vector<int> rso;
    
    std::vector<int> reswapStartOrder;
    int swapPos;
    
    
public:
    
    LinkedList() {
        size.resize(0);
        head.resize(0);
        tail.resize(0);
    }
    
    LinkedList(std::string infile)
    {
        size.resize(0);
        head.resize(0);
        tail.resize(0);
        
        std::ifstream in(infile.c_str());
        std::string s;
        size_t numbreaks = 0;
        size_t numline = 0;
        std::vector<std::vector<std::vector<double> > > dataset;
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
                        std::vector<double> data(n);
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
                        std::vector<std::vector<double> > dset2 = dataset[numline];
                        std::vector<double> nodedat= dset2[std::atoi(ent[1].c_str())];
                        pushBack(nodedat, numline, std::atoi(ent[0].c_str()),std::atoi(ent[1].c_str()));
                    }
                }
                token = s.substr(1,s.length()-2);
                if(numbreaks == 0){
                    size_t n = std::count(token.begin(),token.end(),',')+1;
                    std::vector<double> data(n);
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
                    std::vector<std::vector<double> > dset2 = dataset[numline];
                    std::vector<double> nodedat= dset2[std::atoi(ent[1].c_str())];
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
    
    ~LinkedList() {
        makeCircular(false);
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
    
    void pushBack(T t, int index = 0, int rownum = -1, int colnum = -1) {
        node_ptr linkNode(new LinkNode(t));
        if(rownum == -1)
            linkNode->rownum = index;
        else
            linkNode->rownum = rownum;
        
        if(size.size() <= index){
            size.resize(index+1);
            head.resize(index+1);
            tail.resize(index+1);
        }
        
        ++size[index];
        if(colnum == -1)
            linkNode->colnum = size[index]-1;
        else
            linkNode->colnum = colnum;
        
        if (size[index] == 1) {
            head[index] = tail[index] = linkNode;
            return;
        }
        linkNode->leftNode = tail[index];
        tail[index]->rightNode = linkNode;
        linkNode->leftNode = tail[index];
        tail[index] = linkNode;
    }
    
    void shiftAll(int row, T shift){
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
    
    void setOld(){
        for(int i = 0; i < size.size(); i++){
            node_ptr temp = head[i];
            for(int j = 0; j<size[i];j++){
                temp->olddat = temp->data;
                temp = temp->rightNode;
            }
        }
    }
    
    void resetOld(){
        for(int i = 0; i < size.size(); i++){
            node_ptr temp = head[i];
            for(int j = 0; j<size[i];j++){
                temp->data = temp->olddat;
                temp = temp->rightNode;
            }
        }
    }
    
    void setswap(std::vector<int> i, std::vector<int> j, int pos, int dist = 0){
        swapPos = (pos+dist);
        swapStartOrders.resize(i.size());
        swapEndOrder.resize(j.size());
        
        std::copy(i.begin(), i.end(), swapStartOrders.begin());
        std::copy(j.begin(), j.end(), swapEndOrder.begin());
        
        if(swapPos > size[0])
            reswapStartOrder = glc_per(i);
        else
            reswapStartOrder = i;
    }
    
    void swap(bool reverse = false){
        std::vector<int> i;
        std::vector<int> j;
        i.resize(swapStartOrders.size());
        j.resize(swapEndOrder.size());
        std::copy(swapStartOrders.begin(), swapStartOrders.end(), i.begin());
        std::copy(swapEndOrder.begin(), swapEndOrder.end(), j.begin());
        
        if(!reverse){
            
            if(!i.empty()) {
                int tempPos = 1;
                
                
                std::vector<node_ptr> temps;
                
                for(int ptc = 0; ptc < size.size(); ptc ++)
                    temps.push_back(head[ptc]);
                
                while(swapPos != 0 && tempPos != swapPos){
                    for(std::vector<int>::iterator it = i.begin(); it != i.end(); it++)
                        temps[*it] = temps[*it]->rightNode;
                    ++tempPos;
                }
                
                std::vector<node_ptr> temps2(j.size());
                for(typename std::vector<node_ptr>::iterator it = temps2.begin(); it != temps2.end(); it++)
                    *it = temps[j[it-temps2.begin()]]->rightNode;
                
                for(typename std::vector<int>::iterator it = i.begin(); it != i.end(); it++){
                    if(swapPos != 0){
                        temps[*it]->oldRightNode = temps[*it]->rightNode;
                        temps[*it]->rightNode = temps2[it-i.begin()];
                        temps[*it]->rightNode->leftNode = temps[*it];
                    }
                }
            }
        }
        else{
            if(!i.empty()) {
                int pos = swapPos;
                if(swapPos > size[0]){
                    i = reswapStartOrder;
                    pos = pos%size[0];
                }
                int tempPos = 1;
                
                std::vector<node_ptr> temps;
                
                for(int ptc = 0; ptc < size.size(); ptc ++)
                    temps.push_back(head[ptc]);
                
                while(pos != 0 && tempPos != pos){
                    for(std::vector<int>::iterator it = i.begin(); it != i.end(); it++)
                        temps[*it] = temps[*it]->rightNode;
                    ++tempPos;
                }
                
                for(typename std::vector<int>::iterator it = i.begin(); it != i.end(); it++){
                    if(swapPos != 0){
                        temps[*it]->rightNode = temps[*it]->oldRightNode;
                        temps[*it]->rightNode->leftNode = temps[*it];
                    }
                }
            }
        }
    }
    
    std::vector<int> glc_per(std::vector<int> lc){
        int tempPos = 1;
        
        std::vector<node_ptr> temps;
        
        for(int ptc = 0; ptc < size.size(); ptc ++)
            temps.push_back(head[ptc]);
        
        while(tempPos != size[0]){
            for(std::vector<int>::iterator it = lc.begin(); it != lc.end(); it++)
                temps[*it] = temps[*it]->rightNode;
            ++tempPos;
        }
        
        std::vector<int> nlc;
        for(std::vector<int>::iterator it = lc.begin(); it != lc.end(); it++){
            nlc.push_back(temps[*it]->rightNode->rownum);
        }
        
        return nlc;
    }
    
    std::vector<double> getOne(int row, int slice){
        int tempPos = 0;
        node_ptr temp = head[row];
        
        if((slice < size.size()|| circular) && slice>=0 ){
            while(tempPos != slice){
                temp = temp->rightNode;
                ++tempPos;
            }
            return T(temp->data);
        }
        else if (slice < 0 && circular){
            while(tempPos != slice){
                temp = temp->leftNode;
                --tempPos;
            }
            return T(temp->data);
        }
        else
            return T(0);
    }
    
    void setOne(int row, int slice, T data){
        int tempPos = 0;
        node_ptr temp = head[row];
        if(slice < size.size() || circular){
            while(tempPos != slice){
                temp = temp->rightNode;
                ++tempPos;
            }
            temp->data = data;
        }
    }
    
    std::vector<T> getPair(int row, int start, int dist){
        node_ptr temp;
        node_ptr temp2;
        int tempPos = 0;
        temp = head[row];
        
        while(tempPos != start){
            temp = temp->rightNode;
            ++tempPos;
        }
        temp2 = temp;
        tempPos = 0;
        while(tempPos != dist){
            temp2 = temp2->rightNode;
            ++tempPos;
        }
        std::vector<T> ret(0);
        ret.push_back(temp->data);
        ret.push_back(temp2->data);
        return ret;
    }
    
    std::vector<T> getPairSS(int row1, int row2, int slice){
        node_ptr temp;
        node_ptr temp2;
        int tempPos = 0;
        temp = head[row1];
        temp2 = head[row2];
        
        while(tempPos != slice){
            temp = temp->rightNode;
            temp2 = temp2->rightNode;
            ++tempPos;
        }
        
        std::vector<T> ret(0);
        ret.push_back(temp->data);
        ret.push_back(temp2->data);
        return ret;
        
    }
    
    std::vector<int> getCycles(){
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
    
    
    
    void makeCircular(bool circ = true){
        
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
    
    
    
    void printList(int index = 0){
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
        std::cout << std::endl;
    }
    
    void printListIndices(int index = 0){
        node_ptr temp = head[index];
        if(circular){
            bool headprint = false;
            int numPrint = 0;
            while (!headprint || numPrint < size[index]) {
                if(temp == head[index])
                    headprint = true;
                std::cout << "("<<temp->rownum<<", "<<temp->colnum<<")";
                temp = temp->rightNode;
                numPrint++;
                if(numPrint < size[index])
                    std::cout<< ", ";
            }
            
        }
        else{
            while (temp != NULL) {
                std::cout << "("<<temp->rownum<<", "<<temp->colnum<<"), " << ", ";
                temp = temp->rightNode;
            }
        }
        std::cout << std::endl;
    }
    
    void printListNextIndices(int index = 0){
        node_ptr temp = head[index];
        if(circular){
            bool headprint = false;
            int numPrint = 0;
            while (!headprint || numPrint < size[index]) {
                if(temp == head[index])
                    headprint = true;
                std::cout << "("<<temp->rightNode->rownum<<", "<<temp->rightNode->colnum<<")";
                temp = temp->rightNode;
                numPrint++;
                if(numPrint < size[index])
                    std::cout<< ", ";
            }
            
        }
        else{
            while (temp != NULL) {
                std::cout << "("<<temp->rightNode->rownum<<", "<<temp->rightNode->colnum<<"), " << ", ";
                temp = temp->rightNode;
            }
        }
        std::cout << std::endl;
    }
    
    void printListToFile(int step){
        int numrows = size.size();
        
        std::stringstream sstm;
        sstm << "config_step_" << step <<".txt";
        std::string result = sstm.str();

        std::ofstream out(result.c_str());
        std::streambuf *coutbuf = std::cout.rdbuf();
        std::cout.rdbuf(out.rdbuf());
        for(int j = 0; j < numrows; j++)
            printList(j);
        std::cout << "\n";
        for(int j = 0; j < numrows; j++)
            printListIndices(j);
        std::cout << "\n";
        for(int j = 0; j < numrows; j++)
            printListNextIndices(j);
        std::cout.rdbuf(coutbuf);
    }
    
    
    std::vector<int> getRSO(){
        return rso;
    }
    
    void setPrevSwap(){
        pso = swapStartOrders;
        peo = swapEndOrder;
        rso = reswapStartOrder;
    }
    
    void resetSwap(){
        swapStartOrders = pso;
        swapEndOrder = peo;
        reswapStartOrder = rso;
    }
};


#endif
