/* 
 * File:   adj_handler.h
 * Author: fake
 *
 * Created on May 27, 2012, 8:31 PM
 */

#ifndef ADJ_HANDLER_H
#define	ADJ_HANDLER_H

#include "strutture.h"

#define LEAST -2147483648

/*
 * La matrice di adiacenza e' in realta' un vettore:
 * supponiamo che i primi due siti abbiano vicini 
 * v=[v1,v2,v3]
 * w=[w1,w2,w3,w4]
 * il tutto e' organizzato cosi:
 * adj=[-v1,v2,v3,-w1,w2,w3,w4,-1]
 * gli elementi negativi marcano l'inizio di un altro sito
 */

class neigh_factory{
private:
    int n;  
    int N;
    int buffer[1000];
    int *vicinato[2];
    const int *configuration;
    const int *adj;
    // indice del sito corrente
    int _site;
    // indice degli elementi del vettore di adiacenza corrispondenti a site
    int adj_counter;

    void (neigh_factory::*fetch)(void);
    void f1();
    void f2();
    
public:    
    int operator[](int i) const {
        return buffer[i];
    }
    
    int next(){
        (this->*fetch)();
        return (this->_site);
    }
    
    int site(){
        return this->_site;
    }
    
    int length(){
        return this->n;
    }
    
    void init(const general_partition &p1, const general_partition &p2);
    
    void init(const int *valori_siti, const int *adj1, int N1);
};

#endif	/* ADJ_HANDLER_H */

