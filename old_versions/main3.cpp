/* 
 * Author: Dawid Crivelli
 *
 * Started on January 11, 2012, 2:30 PM
 * 
 * 
 * Version: 3.0, 2012/01/14
 * -Generic partition building
 * -Partitioning of arbitrary symbol strings
 * -Generic intersection and (unreduced) distance
 * -Inheritance classes
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

#include "strutture.h"



options opts;

void general_partition::print() {
    if (opts.verbose > 3) {
        fprintf(stderr, "Adv: {%d", 1);
        for (int i = 1; i < N; i++)
            fprintf(stderr, ",%d", labels[i] - labels[0] + 1);
        fprintf(stderr, "}\n");
    }
    fprintf(stderr, "Sizes: {%d", atom_sizes[0]);
    for (int i = 1; i < n; i++)
        fprintf(stderr, ",%d", atom_sizes[i]);
    fprintf(stderr, "}\n");
}

void partition::print() {
    // {1,3,4,...}
    if (opts.verbose > 2) {
        fprintf(stderr, "Simp: {%d", atom_positions[0]);
        for (int k = 1; k < n; k++)
            fprintf(stderr, ",%d", atom_positions[k]);
        fprintf(stderr, "}\n");
    }
    if (opts.verbose > 3) {
        // {1,0,0,1,0,1,...}
        fprintf(stderr, "{%d", binary[0]);
        for (int k = 1; k < N; k++)
            fprintf(stderr, ",%d", binary[k]);
        fprintf(stderr, "}\n");
    }
}

void entry::print() {
    //sequence name and code
    
    fprintf(stderr, "Seq%03d",this->index);
    if(opts.verbose>1)
        fprintf(stderr, " (%s)",this->name);
    if(this->n < 100)
        fprintf(stderr, ": %s\n", this->seq);
    else
        fprintf(stderr, ": %.20s...\n", this->seq);
    
    //print the whole partition
    if(opts.verbose>2)
        Z->print();
    //adding the entropy information
    fprintf(stderr, "Partitions: %d, Shannon %f, Topological %f\n", Z->n,Z->entropia_shannon(),
            Z->entropia_topologica());
    
    if (X && X->n) {
        //print the whole partition
        if (opts.verbose > 2)
            X->print();
        //adding the entropy information
        fprintf(stderr, "Partitions: %d, Shannon %f, Topological %f\n", X->n, X->entropia_shannon(),
                X->entropia_topologica());
    }
    fprintf(stderr,"\n");
    
}


/*
 * 
 */

int main(int argc, char** argv) {
    int i, j;
    
    double *dist_shan;
    double *dist_shan_r;
    double *dist_top;
    double *dist_top_r;
    double *dist_ham;
    double *dist_fuzzy;
    double *dist_fuzzy_t;
    FILE *out;
    dist_tuple res;
    

    set_program_options(opts,argc,argv);
    
    
    //random number initialization
    srand(opts.seed);
    
    //
    //  ALLOCATION OF MEMORY AND INITIALIZATION
    //
    entry *entries = new entry[opts.n_seq];
    create_temp_arrays(opts);
    

    //distance matrix allocation and zeroing
    dist_shan = new double[opts.n_seq * opts.n_seq];
    dist_shan_r = new double[opts.n_seq * opts.n_seq];
    dist_top = new double[opts.n_seq * opts.n_seq];
    dist_top_r = new double[opts.n_seq * opts.n_seq];
    dist_ham = new double[opts.n_seq * opts.n_seq];
    dist_fuzzy = new double[opts.n_seq * opts.n_seq];
    dist_fuzzy_t = new double[opts.n_seq * opts.n_seq];
    for (i = 0; i < opts.n_seq * opts.n_seq; i++) {
        dist_shan[i] = 0;
        dist_shan_r[i] = 0;
        dist_top[i]=0;
        dist_top_r[i]=0;
        dist_ham[i]=0;
        dist_fuzzy[i]=0;
        dist_fuzzy_t[i]=0;
    }

    
    //
    //  CREATION OF RANDOM SEQUENCES (as entries)
    //
    if(opts.random_sequences){
        fill_entries_randomly(opts,entries);
    }

    //
    //  READ SEQUENCES FROM FILE
    //
    if(opts.random_sequences==false){
        fill_entries_from_file(opts,entries);
    }
    
    if(opts.verbose)
        fprintf(stderr,"Entries filled\n");
    //
    //  ANALISYS OF LOADED SEQUENCES
    //
    for (i = 0; i < opts.n_seq; i++) {
        entries[i].make_partition();
        if (opts.verbose >0)
            entries[i].print();
        
    }
    //
    //  DISTANCE MEASUREMENTS
    //
    if(opts.distance==false)
        exit(0);
    
    if(opts.verbose)
        fprintf(stderr,"Calculating distance matrix\n");
    
    for (i = 0; i < opts.n_seq; i++)
        for (j = i + 1; j < opts.n_seq; j++) {
            int index=i*opts.n_seq+j;
            res = distance(entries[i], entries[j]);
            dist_shan[index] = res.dist_s;
            dist_shan_r[index] = res.dist_s_r;
            dist_top[index]=res.dist_top;
            dist_top_r[index]=res.dist_top_r;
            dist_ham[index]=res.dist_ham;
            
            res = general_distance(entries[i].X,entries[j].X);
            dist_fuzzy[index]=res.dist_fuzzy;
            dist_fuzzy_t[index]=res.dist_fuzzy_t;
        }
    


    //
    //  WRITING THE OUTPUT
    //
    if (opts.write == true) {
        out = fopen("output-dist.bin", "w");
        i = fwrite(dist_shan, sizeof (double), opts.n_seq * opts.n_seq, out);
        fclose(out);

        out = fopen("output-distr.bin", "w");
        i = fwrite(dist_shan_r, sizeof (double), opts.n_seq * opts.n_seq, out);
        fclose(out);

        out = fopen("output-top.bin", "w");
        i = fwrite(dist_top, sizeof (double), opts.n_seq * opts.n_seq, out);
        fclose(out);

        out = fopen("output-topr.bin", "w");
        i = fwrite(dist_top_r, sizeof (double), opts.n_seq * opts.n_seq, out);
        fclose(out);
        
        out = fopen("output-ham.bin", "w");
        i = fwrite(dist_ham, sizeof (double), opts.n_seq * opts.n_seq, out);
        fclose(out);
        
        out = fopen("output-fuzzy.bin", "w");
        i = fwrite(dist_fuzzy, sizeof (double), opts.n_seq * opts.n_seq, out);
        fclose(out);
        
        out = fopen("output-fuzzyt.bin", "w");
        i = fwrite(dist_fuzzy_t, sizeof (double), opts.n_seq * opts.n_seq, out);
        fclose(out);
     
        if(opts.verbose)
            fprintf(stderr,"Written 7x distance matrix\n");
    }
    //
    //  PROGRAM EXIT
    //
    return 0;
}