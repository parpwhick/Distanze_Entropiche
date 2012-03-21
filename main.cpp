/* 
 * Author: Dawid Crivelli
 *
 * Started on January 11, 2012, 2:30 PM
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
#include <ctime>
#include <vector>
#include <map>

#include "strutture.h"


#ifdef _OPENMP
#include <omp.h>
#endif

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
    fprintf(stderr,"media %.1f con %.1f siti/atomo\n",mean,opts.seq_len/mean);
    
}
template void print_partition_stats(linear_partition *, const char* );
template void print_partition_stats(general_partition *, const char* );


void mutation_entropy(std::string *entries){
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
}


int test();

int main(int argc, char** argv) {
    int i=0;
    
    //test();
    
    double *dist_shan;
    double *dist_shan_r;
    double *dist_top;
    double *dist_top_r;
    double *dist_ham;
    double *dist_fuzzy;
    double *dist_fuzzy_t;
    double *dist_fuzzy_r;
    double *dist_fuzzy_r_t;
    FILE *out=0;
        
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
    
    if(opts.from == (FROM_FILE | SEQUENCE)){
        fill_seq_from_file(opts, char_entries);
     //   mutation_entropy(char_entries);
        for (i = 0; i < opts.n_seq; i++){            
        X[i].fill(char_entries[i].data(), opts.seq_len);
        Z[i].from_linear_sequence(char_entries[i].data(), opts.seq_len);
        }        
    }    else if(opts.from == (FROM_FILE | LATTICE)){
        load_lattices_from_file(opts, num_entries);
    
        for (i = 0; i < opts.n_seq; i++)
                Z[i].from_square_lattice(num_entries[i], opts.lato, 2);
    }
    
    if(opts.from & RANDOM){
        for (i = 0; i < opts.n_seq; i++){
          generate_next_sequence(char_entries[0]);
          if(opts.from & LATTICE){
                Z[i].from_square_lattice(char_entries[0].data(), opts.lato, 2);
          }
          else if(opts.from & SEQUENCE){
                X[i].fill(char_entries[0].data(), opts.seq_len);
                Z[i].from_linear_sequence(char_entries[0].data(), opts.seq_len);
          }
        }
    }
       
    printf("Loaded %d sequences long %d\n",opts.n_seq,opts.seq_len);
    if(da_calcolare & SHAN)
        print_partition_stats(X,"semplici");
    if(da_calcolare & GENERAL)
        print_partition_stats(Z,"con salto");
    printf("\n");
    
    
    //
    //  DISTANCE MEASUREMENTS
    //
    if(opts.distance==false)
        exit(0);
    
    
    //distance matrix allocation and zeroing,if needed
    if(da_calcolare & SHAN) dist_shan = new double[opts.n_seq * opts.n_seq];
    if(da_calcolare & RID) dist_shan_r = new double[opts.n_seq * opts.n_seq];
    if(da_calcolare & SHAN_TOP) dist_top = new double[opts.n_seq * opts.n_seq];
    if(da_calcolare & RID_TOP) dist_top_r = new double[opts.n_seq * opts.n_seq];
    if(da_calcolare & HAMM) dist_ham = new double[opts.n_seq * opts.n_seq];
    if(da_calcolare & GENERAL) dist_fuzzy = new double[opts.n_seq * opts.n_seq];
    if(da_calcolare & GENERAL_TOP) dist_fuzzy_t = new double[opts.n_seq * opts.n_seq];
    if(da_calcolare & GENERAL_RID) dist_fuzzy_r = new double[opts.n_seq * opts.n_seq];
    if(da_calcolare & GENERAL_RID_TOP) dist_fuzzy_r_t = new double[opts.n_seq * opts.n_seq];
    
    for (i = 0; i < opts.n_seq * opts.n_seq; i++) {
        if (da_calcolare & SHAN) dist_shan[i]=0;
        if (da_calcolare & RID) dist_shan_r[i]=0;
        if (da_calcolare & SHAN_TOP) dist_top[i]=0;
        if (da_calcolare & RID_TOP) dist_top_r[i]=0;
        if (da_calcolare & HAMM) dist_ham[i]=0;
        if (da_calcolare & GENERAL) dist_fuzzy[i]=0;
        if (da_calcolare & GENERAL_TOP) dist_fuzzy_t[i]=0;
        if(da_calcolare & GENERAL_RID) dist_fuzzy_r[i]=0;
        if(da_calcolare & GENERAL_RID_TOP) dist_fuzzy_r_t[i]=0;

    }

    fprintf(stderr,"Calculating distance matrix\n");
    
    distance d(opts.seq_len);
    
    
    std::clock_t start=std::clock();
    double time_diff;
    double completed_ratio;

#pragma omp parallel for firstprivate(d) schedule(dynamic) num_threads(opts.threads)
    for (i = 0; i < opts.n_seq; i++) {
        for (int j = i + 1; j < opts.n_seq; j++) {
            int index = i * opts.n_seq + j;

            
            if (da_calcolare & (SHAN | RID)) {
                d.binary_partition(X[i], X[j], da_calcolare & RID);
                if (da_calcolare & SHAN) dist_shan[index] = d.dist_s;
                if (da_calcolare & SHAN_TOP) dist_top[index] = d.dist_top;
                if (da_calcolare & RID) dist_shan_r[index] = d.dist_s_r;
                if (da_calcolare & RID_TOP) dist_top_r[index] = d.dist_top_r;
                
//                printf("(%d, %d) semplice: %.3f norm e %.3f ridotta\n",i,j,d.dist_s, d.dist_s_r);
            }

            if (da_calcolare & (GENERAL | GENERAL_RID)) {

                d.fill(Z[i], Z[j]);
                if (da_calcolare & GENERAL) dist_fuzzy[index] = d.dist_fuzzy;
                if (da_calcolare & GENERAL_TOP) dist_fuzzy_t[index] = d.dist_fuzzy_t;
                if (da_calcolare & GENERAL_RID) dist_fuzzy_r[index] = d.dist_fuzzy_r;
                if (da_calcolare & GENERAL_RID_TOP) dist_fuzzy_r_t[index] = d.dist_fuzzy_r_t;
//                printf("(%d, %d) generale: %.3f norm e %.3f ridotta\n\n",i,j,d.dist_fuzzy, d.dist_fuzzy_r);
            }

            if (da_calcolare & HAMM) {

                d.hamming_distance(char_entries[i].c_str(), char_entries[j].c_str());
                dist_ham[index] = d.dist_ham;
            }
           
        }

#ifdef _OPENMP
        int this_thread = omp_get_thread_num();
        if (this_thread)
            continue;
        double time_ratio= omp_get_num_threads();
#else
        double time_ratio = 1.0;
#endif
        fprintf(stderr, "\r");
        time_diff = (std::clock() - start) / (double) CLOCKS_PER_SEC / time_ratio;
		int k=i+1;
        completed_ratio = (2 * k * opts.n_seq - k * k - k + 0.0) / (opts.n_seq * (opts.n_seq - 1));
        fprintf(stderr, "%.1f%% done, ETA %.0fs    ",
                completed_ratio * 100, ceil(time_diff * (1 / completed_ratio - 1)));
        fflush(stderr);

    }
    time_diff = (std::clock() - start) / (double) CLOCKS_PER_SEC;
    fprintf(stderr, "\r100%% done in %.1f seconds of CPU time\n", time_diff);

#define WRITE(where,what) {out=fopen(where, "w");i = fwrite(what, sizeof (double), opts.n_seq * opts.n_seq, out);fclose(out);count++;}
    //
    //  WRITING THE OUTPUT
    //
    if (opts.write == true) {
        int count=0;
        if (da_calcolare & SHAN) 
            WRITE("output-distn.bin",dist_shan);
        
        if (da_calcolare & RID) 
            WRITE("output-distr.bin",dist_shan_r);
        
        if (da_calcolare & SHAN_TOP) 
            WRITE("output-distt.bin",dist_top);
        
        if (da_calcolare & RID_TOP) 
            WRITE("output-distrt.bin",dist_top_r);
        
        if (da_calcolare & GENERAL_RID) 
            WRITE("output-fuzzyr.bin", dist_fuzzy_r);
        
        if (da_calcolare & GENERAL_RID_TOP)
            WRITE("output-fuzzyrt.bin", dist_fuzzy_r_t);
        
        if (da_calcolare & GENERAL) 
            WRITE("output-fuzzy.bin",dist_fuzzy);
        
        if (da_calcolare & GENERAL_TOP) 
            WRITE("output-fuzzyt.bin",dist_fuzzy_t);
        
        if (opts.verbose)
            fprintf(stderr, "Written %dx distance matrix\n",count);
    }
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
