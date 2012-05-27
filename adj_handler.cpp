/* 
 * File:   adj_handler.cpp
 * Author: fake
 * 
 * Created on May 27, 2012, 8:31 PM
 */

#include "adj_handler.h"

extern options opts;

#define nnu(i) (i - (i % lato)+ ((i+lato-1)%lato))
#define nnd(i) ((i/lato)*lato + ((i+lato+1)%lato))
#define nnl(i) (i+N-lato)%N
#define nnr(i) (i+N+lato)%N

typedef struct {
    int *adj;
    int *index;
} adj_struct;

int *adiacenza_square_lattice(int lato){    
    int N=lato*lato;

    int *adj=new int[2*N+1];
    for (int i=0; i<N; i++){
        adj[2*i]   = -nnu(i);
        adj[2*i+1] =  nnl(i);
    }
    adj[2*N]=LEAST;

    return(adj);
}

int *adiacenza_sierpinski(int potenza) {
    int N = 1 << potenza; // lunghezza del lato del triangolo piu grosso = 2^potenza 
    int *line1 = new int[N];
    int *labels1 = new int[N];
    int *line2 = new int[N];
    int *labels2 = new int[N];
    int sitecount = 0;

    int maxsite = 1;
    for (int i = 0; i < potenza; i++)
        maxsite *= 3;

    int *adj = new int[(maxsite - 1)*3+1];
    int adj_count = 0;

    for (int i = 0; i < N; i++)
        line1[i] = line2[i] = 0;
    line1[0] = 1;

    for (int i = 0; i < N; i++) {
        // riempimento dei siti della riga corrente
        line2[0] = 1;
        for (int j = 1; j < i + 1; j++)
            line2[j] = (line1[j] + line1[j - 1]) % 2;

        // labelling dei siti nonnulli
        for (int j = 0; j < i + 1; j++) {
            if (line2[j]) {
                int sign = -1;
                labels2[j] = sitecount;
                sitecount++;
               
                if (i == 0)
                    continue;
                //aggiungi all'elenco dei vicini
                if (j) {
                    if (line1[j - 1]) {
                        adj[adj_count++] = sign * labels1[j - 1];
                        sign = 1;
                    }
                    if (line1[j]) {
                        adj[adj_count++] = sign * labels1[j];
                        sign = 1;
                    }
                    if (line2[j-1]) {
                        adj[adj_count++] = sign * labels2[j-1];
                        sign = 1;
                    }
                } else {
                    if (line1[j])
                        adj[adj_count++] = -labels1[j];
                }

            }
        }
        std::swap(line1, line2);
        std::swap(labels1, labels2);
    }
    adj[(maxsite-1)*3]=LEAST;
    return(adj);
}

int *adiacenza_simple_line(int N){ 
    int *adj=new int[N+1];
    
    adj[0]=LEAST;
    for (int i=1; i<N; i++)
        adj[i]=i-1;
    adj[N]=LEAST;

    return(adj);
}

int *adiacenza_fuzzy_line(int N){
    int *adj=new int[(opts.fuzzy+1)*N];
    int adj_count=0;
    
    adj[0]=LEAST;
    for (int i=1; i<N; i++){
        adj[adj_count++]=-(i-1);
    
        for(int j=2; i-j>=0 && j<=opts.fuzzy+1; j++)
                adj[adj_count++]=i-j;
    }
    adj[N]=LEAST;

    return(adj);
}



void neigh_factory::f1() {
    n = 2;
    buffer[0] = vicinato[0][_site];
    buffer[1] = vicinato[1][_site];
    _site++;
    if (_site >= N)
        _site = -1;
}

void neigh_factory::f2() {
    n = 0;
    int s1;
    int quanti = 1;
    if (adj[adj_counter] == LEAST)
        quanti = 0;

    while (adj[adj_counter + quanti + 1] > 0)
        quanti++;
    for (int i = 0; i < quanti; i++) {
        s1 = adj[adj_counter++];
        if (i == 0)
            s1 = -s1;
        //if(s1 < 0 || s1==site || 
        if (configuration[_site] != configuration[s1])
            continue;
        buffer[n++] = s1;
    }
    _site++;
    if (_site >= N)
        _site = -1;
}

void neigh_factory::init(const general_partition &p1, const general_partition &p2) {
    vicinato[0] = p1.prev_site;
    vicinato[1] = p2.prev_site;
    fetch = &neigh_factory::f1;
    this->N = p1.N;
    _site = 0;
}

void neigh_factory::init(const int *valori_siti, const int *adj1, int N1) {
    configuration = valori_siti;
    fetch = &neigh_factory::f2;
    adj_counter = 0;
    this->N = N1;
    adj = adj1;
    _site = 0;
}


