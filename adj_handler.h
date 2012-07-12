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

/* Dato il link k-esimo, che collega i <--k--> j, 
 * gli array adi[k] e adj[k] restituiscono i e j.
 * 
 * index[i] restituisce il primo link che coinvolge i, fino a index[i+1] sono
 * tutti i link del sito i-esimo, ordinati.
 * 
 * z=fetch(i) restituisce il nr. di coordinazione di 'i' e carica i link opportuni
 * nell'array vicini[0 .. z-1]
 */

class adj_struct {
public:
    int *adi;
    int *adj;
    int *index;
    int N;
    int n_link;
    int zmax;
    mutable int *vicini;
    
    int fetch(int site) const {
        if(site>=N)
            return 0;
        int z=index[site+1]-index[site];
        vicini=adj+index[site];
        __builtin_prefetch(vicini,0,0);
        return(z);
    }
};

adj_struct adiacenza_fuzzy_line(int N);
adj_struct adiacenza_simple_line(int N);
adj_struct adiacenza_square_lattice(int lato);
adj_struct adiacenza_from_file(const char *name_vec1, const char *name_vec2, int & N);
adj_struct adiacenza_sierpinski(int GEN, int &total_size);
#endif	/* ADJ_HANDLER_H */

