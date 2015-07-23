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


class LinkedList {
private:
    
    class LinkNode {
    public:
        typedef std::tr1::shared_ptr<LinkNode> node_ptr;
        
        int rownum;
        T data;
        T olddat;
        node_ptr leftNode;
        node_ptr rightNode;
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
    std::vector<int> reswapStartOrders;
    std::vector<int> reswapEndOrder;
    int swapPos;
    
    
public:
    
    LinkedList() {
        size.resize(0);
        head.resize(0);
        tail.resize(0);
    }
    
    LinkedList(const LinkedList &ll)
    {
        node_ptr temp = node_ptr();
        std::vector<int> oldsize(ll.size);
        for(std::vector<int>::iterator it = oldsize.begin(); it != oldsize.end(); it++){
            temp = ll.head[it-oldsize.begin()];
            for(int i = 0; i < *it; i++){
                pushBack(temp->data,it-oldsize.begin());
                temp = temp->rightNode;
            }
        }
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
    
    void pushBack(T t, int index = 0) {
        node_ptr linkNode(new LinkNode(t));
        linkNode->rownum = index;
        
        if(size.size() <= index){
            size.resize(index+1);
            head.resize(index+1);
            tail.resize(index+1);
        }
        
        ++size[index];
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
        swapPos = (pos+dist)%size[0];
        swapStartOrders.resize(i.size());
        swapEndOrder.resize(j.size());
        
        std::copy(i.begin(), i.end(), swapStartOrders.begin());
        std::copy(j.begin(), j.end(), swapEndOrder.begin());
        
        
        if(!i.empty()){
            
            int max = *std::max_element(i.begin(), i.end());
            
            std::vector<int> identity(max+1);
            iota(identity.begin(),identity.end(),0);
            
            std::vector<int> invj(max+1);
            iota(invj.begin(),invj.end(),0);
            
            for(std::vector<int>::iterator it = j.begin(); it != j.end(); it++)
                invj[*it] = (i[it-j.begin()]);
            
            reswapStartOrders.resize(identity.size());
            reswapEndOrder.resize(identity.size());
            
            std::copy(identity.begin(), identity.end(), reswapStartOrders.begin());
            std::copy(invj.begin(), invj.end(), reswapEndOrder.begin());
        }
        else{
            reswapStartOrders.resize(0);
            reswapEndOrder.resize(0);
        }
        
    }
    void swap(bool reverse = false){
        std::vector<int> i;
        std::vector<int> j;
        
        if(reverse){
            i.resize(reswapStartOrders.size());
            j.resize(reswapEndOrder.size());
            std::copy(reswapStartOrders.begin(), reswapStartOrders.end(), i.begin());
            std::copy(reswapEndOrder.begin(), reswapEndOrder.end(), j.begin());
            
        }
        else{
            i.resize(swapStartOrders.size());
            j.resize(swapEndOrder.size());
            std::copy(swapStartOrders.begin(), swapStartOrders.end(), i.begin());
            std::copy(swapEndOrder.begin(), swapEndOrder.end(), j.begin());
        }
        
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
                    temps[*it]->rightNode = temps2[it-i.begin()];
                    temps[*it]->rightNode->leftNode = temps[*it];
                }
                else if(circular){
                    temps[*it] = temps2[it-i.begin()]->leftNode;
                    temps[*it]->rightNode = temps2[it-i.begin()];
                }
            }
        }
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
    
    
    
    void printLinkedList(int index = 0){
        node_ptr temp = head[index];
        if(circular){
            bool headprint = false;
            int numPrint = 0;
            while (!headprint || numPrint < size[index]) {
                if(temp == head[index])
                    headprint = true;
                T beep = temp->data;
                for(typename T::iterator it = beep.begin(); it != beep.end(); it++)
                    std::cout << *it << "\t";
                temp = temp->rightNode;
                numPrint++;
            }
            
        }
        else{
            while (temp != NULL) {
                T beep = temp->data;
                for(typename T::iterator it = beep.begin(); it != beep.end(); it++)
                    std::cout << *it << "\t";
                temp = temp->rightNode;
            }
        }
        std::cout << std::endl;
    }
    
    void printLinkedListIndices(int index = 0){
        node_ptr temp = head[index];
        if(circular){
            bool headprint = false;
            int numPrint = 0;
            while (!headprint || numPrint < size[index]) {
                if(temp == head[index])
                    headprint = true;
                if(numPrint % size[index] == 0)
                    std::cout << "||\t";
                std::cout << temp->rownum << "\t";
                temp = temp->rightNode;
                numPrint++;
            }
            
        }
        else{
            while (temp != NULL) {
                std::cout << temp->rownum << "\t";
                temp = temp->rightNode;
            }
        }
        std::cout << std::endl;
    }
    
    
    
    void printLinkedList(int index, std::ofstream &f){
        node_ptr temp = head[index];
        if(circular){
            bool headprint = false;
            int numPrint = 0;
            while (!headprint || numPrint < size[index]+1) {
                if(temp == head[index])
                    headprint = true;
                T beep = temp->data;
                for(typename T::iterator it = beep.begin(); it != beep.end(); it++){
                    f << *it << ", ";
                }
                f << std::endl;
                temp = temp->rightNode;
                numPrint++;
            }
            
        }
        else{
            while (temp != NULL) {
                T beep = temp->data;
                for(typename T::iterator it = beep.begin(); it != beep.end(); it++)
                    f << *it << ", ";
                f<<std::endl;
                temp = temp->rightNode;
            }
        }
        f.close();
    }
    
    T returnLinkedList(int index){
        node_ptr temp = head[index];
        T p;
        if(circular){
            bool headprint = false;
            int numPrint = 0;
            while (!headprint || numPrint < size[index]*2) {
                if(temp == head[index])
                    headprint = true;
                T data = temp->data;
                p.push_back(data[0]);
                temp = temp->rightNode;
                numPrint++;
            }
        }
        else{
            while (temp != NULL) {
                T data = temp->data;
                p.push_back(data[0]);
                temp = temp->rightNode;
            }
        }
        return p;
    }
    
};

#endif
