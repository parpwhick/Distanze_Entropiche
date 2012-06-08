/* 
 * File:   adj_handler.h
 * Author: fake
 *
 * Created on May 27, 2012, 8:31 PM
 */

#ifndef ADJ_HANDLER_H
#define	ADJ_HANDLER_H

#include <cstdio>
#include <cstdlib>

#define LEAST -1

/*
 * La matrice di adiacenza e' in realta' un vettore:
 * supponiamo che i primi due siti abbiano vicini 
 * v=[v1,v2,v3]
 * w=[w1,w2,w3,w4]
 * il tutto e' organizzato cosi:
 * adj=[-v1,v2,v3,-w1,w2,w3,w4,-1]
 * gli elementi negativi marcano l'inizio di un altro sito
 */

class adj_struct {
public:
    int *adj;
    int *index;
    int *adi;
    int N;
    int n_link;
    int n;
    int *vicini;
    int zmax;
    
    int fetch(int site){
        if(site>=N)
            return 0;
        n=index[site+1]-index[site];
        vicini=adj+index[site];
        return(n);
    }
};

void ising_simulation(adj_struct NN);
adj_struct adiacenza_fuzzy_line(int N);
adj_struct adiacenza_simple_line(int N);
adj_struct adiacenza_square_lattice(int lato);
adj_struct adiacenza_from_file(const char *name_vec1, const char *name_vec2, int & N);
#endif	/* ADJ_HANDLER_H */

