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
    
    void swap(std::vector<int> i, std::vector<int> j, int pos, int dist = 0){
        
        if(!i.empty()) {
            int tempPos = 1;
            
            std::vector<node_ptr> temps;
            
            for(int ptc = 0; ptc < size.size(); ptc ++)
                temps.push_back(head[ptc]);
            
            while(pos != 0 && tempPos != pos){
                for(std::vector<int>::iterator it = i.begin(); it != i.end(); it++)
                    temps[*it] = temps[*it]->rightNode;
                ++tempPos;
            }
            
            if(dist != 0){
                pos = dist;
                tempPos = 0;
                while(pos != 0 && tempPos != pos){
                    for(std::vector<int>::iterator it = i.begin(); it != i.end(); it++)
                        temps[*it] = temps[*it]->rightNode;
                    ++tempPos;
                }
                
            }
            
            std::vector<node_ptr> temps2(j.size());
            for(typename std::vector<node_ptr>::iterator it = temps2.begin(); it != temps2.end(); it++)
                *it = temps[j[it-temps2.begin()]]->rightNode;
            
            for(typename std::vector<int>::iterator it = i.begin(); it != i.end(); it++){
                if(pos != 0){
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
    
    void reverseswap(std::vector<int> i, std::vector<int> j, int pos, int dist = 0){
        
        if(!i.empty()){
            std::vector<int> identity(i.back()+1);
            iota(identity.begin(),identity.end(),0);
            
            std::vector<int> invj(i.back()+1);
            iota(invj.begin(),invj.end(),0);
            
            for(std::vector<int>::iterator it = j.begin(); it != j.end(); it++)
                invj[*it] = (i[it-j.begin()]);
            
            swap(identity,invj,pos,dist);
        }
        
    }
    
    
    std::vector<double> getOne(int row, int slice){
        int tempPos = 0;
        node_ptr temp = head[row];
        while(tempPos != slice){
            temp = temp->rightNode;
            ++tempPos;
        }
        return T(temp->data);
    }
    
    void setOne(int row, int slice, T data){
        int tempPos = 0;
        node_ptr temp = head[row];
        while(tempPos != slice){
            temp = temp->rightNode;
            ++tempPos;
        }
        temp->data = data;
    }
    
    std::vector<T> getPair(int row, int start, int dist){
        node_ptr temp;
        node_ptr temp2;
        if(dist!=0){
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
        }
        else{
            int tempPos = 0;
            temp = head[row];
            temp2 = head[start];
            while(tempPos != start){
                temp = temp->rightNode;
                temp2 = temp2->rightNode;
                ++tempPos;
            }
        }
        std::vector<T> ret(0);
        ret.push_back(temp->data);
        ret.push_back(temp2->data);
        return ret;        
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
            while (!headprint || numPrint < size[index]+3) {
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
    
    void printLinkedList(int index, std::ofstream &f){
        node_ptr temp = head[index];
        if(circular){
            bool headprint = false;
            int numPrint = 0;
            while (!headprint || numPrint < size[index]*2) {
                if(temp == head[index])
                    headprint = true;
                T beep = temp->data;
                for(typename T::iterator it = beep.begin(); it != beep.end(); it++)
                    f << *it << ", ";
                temp = temp->rightNode;
                numPrint++;
            }
            
        }
        else{
            while (temp != NULL) {
                T beep = temp->data;
                for(typename T::iterator it = beep.begin(); it != beep.end(); it++)
                    f << *it << ", ";
                temp = temp->rightNode;
            }
        }
        f << std::endl;
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
