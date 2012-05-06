/* 
 * Author: Dawid Crivelli
 *
 * Started on January 11, 2012, 2:30 PM
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
#include <cstring>
#include <string>
#include <vector>
#include <map>
#include <time.h>
#include "strutture.h"

extern options opts;
double *mylog=0;
int *colore=0;



template <typename part>
void print_partition_stats(part *X, const char* name){
    int min=X[0].N*20, max=0;
    double mean=0, std=0;
    for(int i=0; i<opts.n_seq; i++){
        min=std::min(X[i].n,min);
        max=std::max(X[i].n,max);
        mean+=X[i].n;
        std+=(X[i].n)*(X[i].n);
        
    }
    mean /= opts.n_seq;
    std /= opts.n_seq;
    std = std - mean*mean;
    
    fprintf(stderr,"Partizioni %9s: nr. frammenti tra [%d,%d], ",name,min,max);
    fprintf(stderr,"media %.2f con %.2f siti/atomo\n",mean,opts.seq_len/mean);
    
}
template void print_partition_stats(linear_partition *, const char* );
template void print_partition_stats(general_partition *, const char* );


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
            
    set_program_options(opts,argc,argv);
    
    int &da_calcolare=opts.da_calcolare;
        
    //random number initialization
    srand(time(0));
    xrandinit(time(0));    
    
    //
    //  ALLOCATION OF MEMORY AND INITIALIZATION
    //
    linear_partition *X = new linear_partition[opts.n_seq];
    general_partition *Z = new general_partition[opts.n_seq];
    
    //logarithm lookup table, 6x program speedup
    int lunghezza=std::max(opts.seq_len, opts.n_seq) + 10;
    mylog = new double[ lunghezza ];
    for (int i = 1; i <  lunghezza;  i++)
        mylog[i] = log(i);
    mylog[0]=0;
    
    
    std::string *char_entries=new std::string[opts.n_seq];
    int **num_entries=new int*[opts.n_seq];

    if (opts.topologia == LINEARE && opts.letto_da == FROM_FILE) {
        fill_seq_from_file(opts, char_entries);
        //   mutation_entropy(char_entries);
        for (int i = 0; i < opts.n_seq; i++) {
            X[i].fill(char_entries[i].data(), opts.seq_len);
            Z[i].from_linear_sequence(char_entries[i].data(), opts.seq_len);
        }
    } else if (opts.topologia == RETICOLO && opts.letto_da == FROM_FILE ) {
        load_lattices_from_file(opts, num_entries);

        for (int i = 0; i < opts.n_seq; i++)
            Z[i].from_square_lattice(num_entries[i], opts.lato, 2);
    }

    if (opts.letto_da== RANDOM) {
        for (int i = 0; i < opts.n_seq; i++) {
            generate_next_sequence(char_entries[0]);
            if (opts.topologia == RETICOLO) {
                Z[i].from_square_lattice(char_entries[0].data(), opts.lato, 2);
            } else if (opts.topologia == LINEARE) {
                X[i].fill(char_entries[0].data(), opts.seq_len);
                Z[i].from_linear_sequence(char_entries[0].data(), opts.seq_len);
            }
        }
    }

    printf("Loaded %d sequences long %d\n", opts.n_seq, opts.seq_len);
    if (da_calcolare & SHAN)
        print_partition_stats(X, "semplici");
    if (da_calcolare & GENERAL)
        print_partition_stats(Z, "con salto");
    printf("\n");


    //
    //  DISTANCE MEASUREMENTS
    //
    if (opts.distance == false)
        exit(0);

    calcola_matrice_distanze(X, Z, char_entries);
    //
    //  PROGRAM EXIT
    //
    delete []X;
    delete []Z;
    delete []mylog;
    delete []char_entries;
    delete []num_entries;
    return 0;
}
