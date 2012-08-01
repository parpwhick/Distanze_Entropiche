/* 
 * Author: Dawid Crivelli
 *
 * Started on January 11, 2012, 2:30 PM
 * 
 * Version: 6.1, 2012/06/21
 * -Added simulation class
 * -Improved entropy calculation, 2x speedup
 * -New random number generator, for simulation use
 * -Cleanups and memory usage reduction
 * 
 * Version: 6.0, 2012/06/01
 * -General adjacency vector input
 * -Removed hashing and unused functions
 *   
 * Version: 5.2, 2012/03/22
 * -Iterators for the linked list
 * -General cleanups
 * 
 * Version: 5.1, 2012/03/11
 * -Super optimized reduction
 * -Hashing implemented everywhere
 * -Reduced number of algorithms
 * -Multidimensional inputs
 * 
 * Version: 5.0, 2012/03/05
 * -Fully working general reduction
 * -with interchangable algorithms
 * 
 * Version: 4.1, 2012/01/27
 * -Common partition factor by percolation
 * -Similarity distance
 *
 * Version: 4.0, 2012/01/20
 * -Multithreading
 * -Algorithm optimization and automatic selection
 *  
 * Version: 3.0, 2012/01/15
 * -Generic partition building
 * -Partitioning of arbitrary symbol strings
 * -Generic intersection and (unreduced) distance
 *  
 * Version: 2.0, 2012/01/13
 * -Reading from files
 * -Help system 
 * -Thorough input options
 * 
 * Version: 1.0, 2012/01/12
 * -Partitioning
 * -Distance matrix calculation
 */

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <map>
#include "strutture.h"
#include "partizioni.h"
#include "distance.h"

extern options opts;
double *mylog=0;

/************************************************
 ************************************************
 *             PARTIZIONI SEMPLICI              *
 ************************************************
 ************************************************
 */

void linear_partition::print() {
    if (opts.verbose > 3) {
        // {1,0,0,1,0,1,...}
        print_array(binary, N, "Binary");
    }
    fprintf(stdout, "Partitions[n]: %d, Shannon %f, Topological %f\n", n, entropia_shannon,
            entropia_topologica);
}

template <typename T>
linear_partition::linear_partition(const T* seq, int len) {
    this->fill(seq, len);
}

template <typename T>
void linear_partition::fill(const T* seq, int len) {
    //Total length of the partition is equal to the sequence
    N = len;
    binary = new int[len];

    //first one always start an atom
    binary[0] = 1;

    //checking for a different symbol from the one before
    //when that happens => new atom!
    for (int i = 1; i < len; i++)
        binary[i] = seq[i] != seq[i - 1];

    auto entropie = ordered_vector_entropy(seq,N);
    n=entropie.second;
    entropia_shannon = entropie.first;
    entropia_topologica = mylog[n];
}
template void linear_partition::fill(const char *, int);

void print_partition_stats(linear_partition *X){
    label_t min=X[0].N*20, max=0;
    double mean=0, std=0;
    for(int i=0; i<opts.n_seq; i++){
	if(X[i].n < 10)
		printf("Few atoms in partition: %d\n",i);
        min=std::min(X[i].n,min);
        max=std::max(X[i].n,max);
        mean+=X[i].n;
        std+=(X[i].n)*(X[i].n);
        
    }
    mean /= opts.n_seq;
    std /= opts.n_seq;
    std = std - mean*mean;
    
    fprintf(stderr,"Partitions: n. atoms between [%d,%d], ",min,max);
    fprintf(stderr,"avg %.2f with %.2f sites/atom\n",mean,opts.seq_len/mean);
    
}

double *mutation_entropy(std::string *entries){
    typedef std::map<char,int> map_t;
    map_t *histogram=new map_t[opts.seq_len];
    
    for(int i=0; i < opts.n_seq; i++){
        for(int j=0; j < opts.seq_len; j++){
            histogram[j][ entries[i][j] ]++;
        }
    }
        
    double *H=new double[opts.seq_len];
    for (int j = 0; j < opts.seq_len; j++) {
        int labels=0;
        H[j]=0;
        for (map_t::iterator ii = histogram[j].begin(); ii != histogram[j].end(); ++ii) {
            H[j] += ii->second * mylog[ii->second];
            labels++;
        }
        H[j]= -H[j]/opts.n_seq+ mylog[opts.n_seq];
    }
//    
//    printf("Entropie sito per sito: {");
//    for (int j = 0; j < opts.seq_len; j++) 
//        printf("%02.1f,",H[j]);
//    printf("}\n");

	return(H);
}

int main(int argc, char** argv) {
            
    opts.partition_type = LINEAR_PARTITION;
    set_program_options(opts,argc,argv);
    
    //
    //  ALLOCATION OF MEMORY AND INITIALIZATION
    //
    linear_partition *X = new linear_partition[opts.n_seq];
    
    //logarithm lookup table, 6x program speedup
    int lunghezza=std::max(opts.seq_len, opts.n_seq) + 10;
    mylog = new double[ lunghezza ];
    for (int i = 1; i <  lunghezza;  i++)
        mylog[i] = log(i);
    mylog[0]=0;
        
    std::string *buffer=new std::string[opts.n_seq];
    if (opts.letto_da == FROM_FILE)
        fill_seq_from_file(opts, buffer);
    else
        fill_entries_randomly(buffer);    
   //   mutation_entropy(char_entries);

    for (int i = 0; i < opts.n_seq; i++) 
        X[i].fill(buffer[i].data(), opts.seq_len);

    printf("Loaded %d sequences long %d\n", opts.n_seq, opts.seq_len);
    print_partition_stats(X);
    printf("\n");

    //
    //  DISTANCE MEASUREMENTS
    //
    if (opts.distance == false)
        exit(0);

    calcola_matrice_distanze(X);
    //
    //  PROGRAM EXIT
    //
    delete []X;
    delete []mylog;
    delete []buffer;
    return 0;
}
