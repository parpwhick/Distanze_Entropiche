#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <stdio.h>
#include <string>
#include <map>
#include <vector>
//#include <unordered_map>

#include "strutture.h"

extern options opts;
extern double *mylog;

void distance::allocate(int n) {
    if (opts.da_calcolare & (SHAN | RID)) {
        common_factor = new int[n];
        reduced1 = new int[n];
        reduced2 = new int[n];
        product_reduced = new int[n];
    }
    product = new u_int64_t[n];
}

distance::distance(int n) {
    N = n;
    allocate(N);
}

distance::distance(const distance &d1) {
    N = d1.N;
    allocate(N);
}

distance::~distance() {
    if (N) {
        if (opts.da_calcolare & (SHAN | RID)) {
            delete[] common_factor;
            delete[] reduced1;
            delete[] reduced2;
            delete[] product_reduced;
        }
        delete[] product;
    }
    N = 0;
}

void distance::fill(const linear_partition& e1, const linear_partition& e2) {
    this->binary_partition(e1, e2);
}

void distance::fill(const general_partition& e1, const general_partition& e2) {
    static int imagecount = 0;

    if (opts.da_calcolare & GENERAL_RID) {
        partizione_comune.linear_intersection(e1, e2);
        //partizione_comune.trivial(e1.N);
        ridotto1.reduce(e1, partizione_comune);
        ridotto2.reduce(e2, partizione_comune);

     
        if(!opts.verbose){
        int quanto=std::min(e1.N,50);
        print_array(partizione_comune.labels,quanto,"lbls comune");
        print_array(e1.labels,quanto,               "lbls e1    ");
        print_array(ridotto1.labels,quanto,         "lbls ridot1");
        print_array(e2.labels,quanto,               "lbls e2    ");   
        print_array(ridotto2.labels,quanto,         "lbls ridot2");
        
        }
//        ridotto1.reduce(e1, e2);
//        ridotto2.reduce(e2, e1);

        linear_product_sorted(ridotto1,ridotto2);
        //linear_product_sorted(e1,e2);
        dist_fuzzy_r = dist_fuzzy;
        dist_fuzzy_r_t = dist_fuzzy_t;
       //printf("comune generale: %d, r1: %d/%d, r2: %d/%d, prod: %d\n",partizione_comune.n, ridotto1.n, e1.n, ridotto2.n,e2.n,(int)dist_fuzzy_t);
        if (opts.graphics && (opts.from & LATTICE)) {
            char filename[255];
            imagecount++;
            sprintf(filename, "ridotto%03d.ppm", imagecount);
            ppmout2(e1.labels, e2.labels, e1.lato, filename);
            imagecount++;
            sprintf(filename, "ridotto%03d.ppm", imagecount);
            ppmout2(ridotto1.labels, ridotto2.labels, e1.lato, filename);
        }        
    }

    if (opts.da_calcolare & GENERAL) {
        switch (opts.alg) {
            case PMATRIX:
            case SORTED:
            default:
                linear_product_sorted(e1, e2);
                break;
        }
    }
}

void print_binary_partition(int*p, int N) {
    // |...|...||...|...
    for (int i = 0; i < N; i++)
        printf("%c", (p[i]) ? '|' : '.');
    printf("\n");
}

template <typename T>
double entropy_binary_partition(const T *p, int n) {
    //p: binary partition
    //n: total length
    int i = 0;
    int mu;
    int begin;
    double H = 0;

    //the first position always starts an atom
    begin = 0;
    for (i = 1; i < n; i++) {
        //whenever we find 1 (a new atom)
        if (p[i]) {
            //the closed (old)atom's length is calculated
            mu = i - begin;
            //the new one is ready to go
            begin = i;
            //we add the entropy, with the trick mu>0 and when mu=1 the log is 0
            if (mu > 1)
                H += (double) mu * mylog[mu];
        }
    }
    //we check the last one, in case it's left hanging
    mu = n - begin;
    if (mu > 1) H += mu * mylog[mu];

    //proper entropy normalization
    H = -H / n + mylog[n];
    return (H);
}
template double entropy_binary_partition(int const*, int);


template <typename T>
void distance::hamming_distance(const T* seq1, const T* seq2) {
    this->dist_ham = 0;
    for (int i = 0; i < N; i++)
        this->dist_ham += (double) (seq1[i] != seq2[i]);
}
template void distance::hamming_distance(const char*,const char*);



//function giving the distance between 2 partitions, by intersecting 
//and calculating the relevant entropy
void distance::binary_partition(const linear_partition &first, const linear_partition &second,int ridotta) {

    int N = first.N;
    
    //the partition product/intersection, OR'ing
    for (int i = 0; i < N; i++)
        product[i] = first.binary[i] | second.binary[i];
    product[0] = 1;
    
    //calculating ALL the entropies (3 per non-reduced Rohlin dist, 3 per reduced)
    double
	    h1 = first.entropia_shannon,
            h2 = second.entropia_shannon,
            h12 = entropy_binary_partition(product, N);

    this->dist_s = 2 * h12 - h1 - h2;
    
    
    //DISTANZA TOPOLOGICA
    int coperture1 = 0, coperture2 = 0, coperture12 = 0;
    for (int i = 0; i < N; i++) {
        coperture1 += first.binary[i];
        coperture2 += second.binary[i];
        coperture12 += product[i];
    }
    this->dist_top = 2 * mylog[coperture12] - mylog[coperture1] - mylog[coperture2];
    
    
    if(!ridotta)
        return;
    //DISTANZA RIDOTTA
    //the common factor of the partition, AND'ing
    for (int i = 0; i < N; i++)
        common_factor[i] = first.binary[i] & second.binary[i];
    common_factor[0] = 1;

    //reducing both the partitions XOR'ing with the common factor
    for (int i = 0; i < N; i++)
        reduced1[i] = common_factor[i] ^ first.binary[i];
    reduced1[0] = 1;
    for (int i = 0; i < N; i++)
        reduced2[i] = common_factor[i] ^ second.binary[i];
    reduced2[0] = 1;

    //intersection of the reduced partitions, OR'ing
    for (int i = 0; i < N; i++)
        product_reduced[i] = reduced1[i] | reduced2[i];
    //  alternatively, one could XOR the unreduced partitions
    //  intersection_reduced[i] =first->binary[i] ^ second->binary[i];
    product_reduced[0] = 1;
    
    double  hr1 = entropy_binary_partition(reduced1, N),
            hr2 = entropy_binary_partition(reduced2, N),
            hr12 = entropy_binary_partition(product_reduced, N);
    this->dist_s_r = 2 * hr12 - hr1 - hr2;
    
    
    //DISTANZA TOPOLOGICA RIDOTTTA
    int coperture1r = 0, coperture2r = 0, coperture12r = 0,coperture_common=0;
    for (int i = 0; i < N; i++) {
        coperture1r += reduced1[i];
        coperture2r += reduced2[i];
        coperture12r += product_reduced[i];
        coperture_common+=common_factor[i];
    }
    this->dist_top_r = 2 * mylog[coperture12r] - mylog[coperture1r] - mylog[coperture2r];
    
    //printf("comune semplice: %d, r1: %d/%d, r2: %d/%d, prod: %d\n",coperture_common, coperture1r, coperture1, coperture2r, coperture2, coperture12r);
}


inline int compare (const void * a, const void * b){
  return ( *(u_int64_t*)a - *(u_int64_t*)b );
}


void distance::linear_product_sorted(const general_partition &p1,const  general_partition &p2) {
    int i;
    int label_count = 0;
    double H = 0;
    int mu;
    int begin;

    for (i = 0; i < N; i++){
        u_int64_t temp1=p1.labels[i];
        u_int64_t temp2=p2.labels[i];
    
        product[i] = (temp1<<32) | temp2;
    }
    
    static int imagecount=0; 
    char filename[255];
    if (opts.graphics && (opts.from & LATTICE)) {
          imagecount++;
          sprintf(filename,"prodotto%03d.ppm",imagecount);
          ppmout2(p1.labels,p2.labels,p1.lato,filename);
          imagecount++;
          sprintf(filename,"prodotto%03d.ppm",imagecount);
          ppmout(product,p1.lato,filename);
      }
    
    qsort(product,N,sizeof(u_int64_t),compare);

    //the first position always starts an atom
    begin = 0;
    label_count=1;
    u_int64_t old_val=product[0];
    for (i = 1; i < N; i++) {
        //whenever we find 1 (a new atom)
        if (product[i]!=old_val) {
            //a new atom starts
            label_count++;
            //the closed (old)atom's length is calculated
            mu = i - begin;
            //the new one is ready to go
            begin = i;
            //cache the new label to check
            old_val=product[i];
            //we add the entropy, with the trick mu>0 and when mu=1 the log is 0
            if (mu > 1)
                H += (double) mu * mylog[mu];
        }
    }
    //we check the last one, in case it's left hanging
    mu = N - begin;
    H += mu * mylog[mu];
    
    //normalize the result
    H = -H / N + mylog[N];

    double h1 = p1.entropia_shannon,
            h2 = p2.entropia_shannon,
            t1 = p1.entropia_topologica,
            t2 = p2.entropia_topologica;
    this->dist_fuzzy = 2 * H - h1 - h2;
    this->dist_fuzzy_t = 2 * mylog[label_count] - t1 - t2;


}