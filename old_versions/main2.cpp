
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>


#include "strutture2.h"


int *common_factor;
int *reduced1;
int *reduced2;
int *product;
int *product_reduced;
double *mylog;


options opts;


double entropy_binary_partition(int *p, int n) {
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

double partition::entropia_shannon() {
    //this is just a helper function
    return ( entropy_binary_partition(this->binary, this->N));
}

double partition::entropia_topologica(){
    return (mylog[this->n]);
}

void entry::make_partition() {
    int i, j;

    //Total length of the partition is equal to the sequence
    Z.N = n;
    Z.binary = new int[n];
    Z.atom_positions = new int[n];

    //first one always start an atom
    Z.binary[0] = 1;

    //checking for a different symbol from the one before
    //when that happens => new atom!
    for (i = 1; i < n; i++)
        Z.binary[i] = seq[i] != seq[i - 1];

    //now going back and writing the atom positions
    j = 0;
    for (i = 0; i < n; i++)
        if (Z.binary[i]) {
            Z.atom_positions[j++] = i + 1;
        }
    //number of atoms found - index++ of last the atom in the array (zero addressing)
    Z.n = j;

}

void print_binary_partition(int*p, int N) {
    // |...|...||...|...
    for (int i = 0; i < N; i++)
        printf("%c", (p[i]) ? '|' : '.');
    printf("\n");
}

//function giving the distance between 2 partitions, by intersecting 
//and calculating the relevant entropy
dist_tuple intersect(partition &first, partition &second) {

    dist_tuple res;
    int N = first.N;
    
    int coperture1=0, coperture2=0,coperture12=0,coperture1r=0,coperture2r=0,coperture12r=0;

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

    //the partition product/intersection, OR'ing
    for (int i = 0; i < N; i++)
        product[i] = first.binary[i] | second.binary[i];
    product[0] = 1;

    //intersection of the reduced partitions, OR'ing
    for (int i = 0; i < N; i++)
        product_reduced[i] = reduced1[i] | reduced2[i];
    //  alternatively, one could XOR the unreduced partitions
    //  intersection_reduced[i] =first.binary[i] ^ second.binary[i];
    product_reduced[0] = 1;

    for (int i = 0; i < N; i++) {
        coperture1 += first.binary[i];
        coperture2 += second.binary[i];
        coperture12+=product[i];
        coperture1r+=reduced1[i];
        coperture2r+=reduced2[i];
        coperture12r += product_reduced[i];
    }
    

    //calculating ALL the entropies (3 per non-reduced Rohlin dist, 3 per reduced)
    double
    h1 = entropy_binary_partition(first.binary, N),
            h2 = entropy_binary_partition(second.binary, N),
            hr1 = entropy_binary_partition(reduced1, N),
            hr2 = entropy_binary_partition(reduced2, N),
            h12 = entropy_binary_partition(product, N),
            hr12 = entropy_binary_partition(product_reduced, N);

    //packing the output tuple with the distances
    res.dist = 2 * h12 - h1 - h2;
    res.dist_r = 2 * hr12 - hr1 - hr2;
    res.dist_top = 2*mylog[coperture12]-mylog[coperture1]-mylog[coperture2] ;
    res.dist_top_r = 2*mylog[coperture12r]-mylog[coperture1r]-mylog[coperture2r];

    //printing the pretty partition graphs if selected
    if (opts.verbose > 2) {
        printf("H(1):\t\t%f\n", h1);
        print_binary_partition(first.binary, N);
        printf("H(2):\t\t%f\n", h2);
        print_binary_partition(second.binary, N);
        if (opts.verbose > 3) {
            printf("Hr(1):\t\t%f\n", hr1);
            print_binary_partition(reduced1, N);
            printf("Hr(2):\t\t%f\n", hr2);
            print_binary_partition(reduced2, N);
        }
        printf("H(12):\t\t%f\n", h12);
        print_binary_partition(product, N);
        printf("Hr(12):\t\t%f\n", hr12);
        print_binary_partition(product_reduced, N);
        printf("\n");

        printf("d(1,2)=\t\t%f\n", res.dist);
        printf("dr(1,2)=\t%f\n", res.dist_r);
        printf("dt(1,2)=\t%f\n", res.dist_top);
        printf("dtr(1,2)=\t%f\n", res.dist_top_r);
    }

    return (res);
}

void partition::print() {
    // {1,3,4,...}
    fprintf(stderr, "{");
    for (int k = 0; k < n; k++)
        fprintf(stderr, "%d,", atom_positions[k]);
    fprintf(stderr, "}\n");

    // {1,0,0,1,0,1,...}
    fprintf(stderr, "{");
    for (int k = 0; k < N; k++)
        fprintf(stderr, "%d,", binary[k]);
    fprintf(stderr, "}\n");

}

void entry::print() {
    //sequence name and code
    if(this->n < 100)
        fprintf(stderr, "Seq%d: %s\n", this->index,this->seq);
    else{
        fprintf(stderr, "Seq%d: %.20s...\n", this->index,this->seq);
    }
    //print the whole partition
    if(opts.verbose>1)
        Z.print();
    //adding the entropy information
    fprintf(stderr, "Partitions: %d, Shannon %f, Topological %f\n\n", Z.n,Z.entropia_shannon(),
            Z.entropia_topologica());
}

/*
 * 
 */

int main(int argc, char** argv) {
    int i, j;
    
    double *distances;
    double *distances_r;
    double *dist_top;
    double *dist_top_r;
    FILE *out;
    dist_tuple res;
    

    set_program_options(opts,argc,argv);
    
    //
    //  ALLOCATION OF MEMORY AN INITIALIZATION
    //
    entry *entries = new entry[opts.n_seq];

    //arrays needed to intersect the partitions
    common_factor = new int[opts.seq_len];
    reduced1 = new int[opts.seq_len];
    reduced2 = new int[opts.seq_len];
    product = new int[opts.seq_len];
    product_reduced = new int[opts.seq_len];

    //logarithm lookup table, 6x program speedup
    mylog = new double[opts.seq_len + 10];
    for (i = 1; i < opts.seq_len + 10; i++)
        mylog[i] = log(i);
    mylog[0]=0;

    //random number initialization
    srand(opts.seed);

    //distance matrix initialization and zeroing
    distances = new double[opts.n_seq * opts.n_seq];
    distances_r = new double[opts.n_seq * opts.n_seq];
    dist_top = new double[opts.n_seq * opts.n_seq];
    dist_top_r = new double[opts.n_seq * opts.n_seq];
    for (i = 0; i < opts.n_seq * opts.n_seq; i++) {
        distances[i] = 0;
        distances_r[i] = 0;
        dist_top[i]=0;
        dist_top_r[i]=0;
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
    
    for (i = 0; i < opts.n_seq; i++)
        for (j = i + 1; j < opts.n_seq; j++) {
            res = intersect(entries[i].Z, entries[j].Z);
            distances[i * opts.n_seq + j] = res.dist;
            distances_r[i * opts.n_seq + j] = res.dist_r;
            dist_top[i*opts.n_seq+j]=res.dist_top;
            dist_top_r[i*opts.n_seq+j]=res.dist_top_r;
        }


    //
    //  WRITING THE OUTPUT
    //
    if (opts.write == true) {
        out = fopen("output-dist.bin", "w");
        i = fwrite(distances, sizeof (double), opts.n_seq * opts.n_seq, out);
        fclose(out);

        out = fopen("output-distr.bin", "w");
        i = fwrite(distances_r, sizeof (double), opts.n_seq * opts.n_seq, out);
        fclose(out);

        out = fopen("output-top.bin", "w");
        i = fwrite(dist_top, sizeof (double), opts.n_seq * opts.n_seq, out);
        fclose(out);

        out = fopen("output-topr.bin", "w");
        i = fwrite(dist_top_r, sizeof (double), opts.n_seq * opts.n_seq, out);
        printf("Written %d records\n", i);
        fclose(out);
    }
    //
    //  PROGRAM EXIT
    //
    return 0;
}


