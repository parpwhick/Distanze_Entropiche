#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>

#include "strutture.h"



int *common_factor;
int *reduced1;
int *reduced2;
int *product;
int *product_reduced;
double *mylog;


extern options opts;

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

template <typename T>
partition::partition(const T* seq, int len) {
    int i, j;
    
    //Total length of the partition is equal to the sequence
    N = len;
        
    binary = new int[len];
    atom_positions = new int[len];

    //first one always start an atom
    binary[0] = 1;

    //checking for a different symbol from the one before
    //when that happens => new atom!
    for (i = 1; i < len; i++)
        binary[i] = seq[i] != seq[i - 1];

    //now going back and writing the atom positions
    j=0;
    for (i = 0; i < len; i++)
        if (binary[i]) {
            atom_positions[j++] = i + 1;
        }
    //number of atoms found - index++ of last the atom in the array (zero addressing)
    n = j;

}

template <typename T>
general_partition::general_partition(const T* seq, int len) {
    int i, j;
    //is not reset between function calls, it starts where it was left the previous time
    static int running_label=0;
    //this run's  atoms count
    int local_count;
    //position of the last site belonging to the atom we're trying to add to
    int last_good;
    
    //Total length of the partition is equal to the sequence
    N = len;
    
    //position of first symbol (0 based index)
    atom_positions = new int[N];
    
    //size of atoms (0 based index)
    atom_sizes=new int[N];
    labels=new int[N];

    //presetting labels for this partition to 0
    for (i=0;i<N;i++)
        labels[i]=0;
    
    //beginning count of identified atoms
    local_count=0;
    //we check every site for belonging to a partition set
    for (i=0;i<N;i++){
        //if the site was already "colored" check the next one
        if (labels[i])
            continue;
        //we give a new label to the site
        ++running_label;
        atom_positions[local_count]=i;
        labels[i]=running_label;
        
        //we know that it belongs (first, trivial)
        last_good=i;
        int this_atom_size=1;
        //now we add all the other possible sites to it, starting from i+1!
        //staying below N
        //and allowing a spacing of at most opts.fuzzy spaces between them
        for(j=i+1; j<last_good+opts.fuzzy+2 && j<N; j++){
            
            //if it belongs to the same cluster...
            if(seq[j]==seq[i]){
                //..we label it accordingly
                labels[j]=running_label;
                //and marking it as "good" (since we're on an ordered line)
                last_good=j;
                this_atom_size++;
            }
        }
        atom_sizes[local_count]=this_atom_size;
        local_count++;
    }
    //this many different atoms were found
    n=local_count;

}

void entry::make_partition(){
    
        Z=new partition(seq,n);
        X=new general_partition(seq,n);
    }


double general_partition::entropia_shannon(){
    double H=0;
    for(int i=0;i<n;i++)
     H += (double) atom_sizes[i] * mylog[atom_sizes[i]];
    H=-H/N+mylog[N];
    return(H);
}



void create_temp_arrays(options opts){
    //arrays needed to intersect the partitions
    common_factor = new int[opts.seq_len];
    reduced1 = new int[opts.seq_len];
    reduced2 = new int[opts.seq_len];
    product = new int[opts.seq_len];
    product_reduced = new int[opts.seq_len];

    //logarithm lookup table, 6x program speedup
    mylog = new double[opts.seq_len + 10];
    for (int i = 1; i < opts.seq_len + 10; i++)
        mylog[i] = log(i);
    mylog[0]=0;
}

void print_binary_partition(int*p, int N) {
    // |...|...||...|...
    for (int i = 0; i < N; i++)
        printf("%c", (p[i]) ? '|' : '.');
    printf("\n");
}

//function giving the distance between 2 partitions, by intersecting 
//and calculating the relevant entropy
dist_tuple distance(entry &e1, entry &e2) {

    partition *first=e1.Z;
    partition *second=e2.Z;
    dist_tuple res;
    int N = e1.n;
    int coperture1=0, coperture2=0,coperture12=0,coperture1r=0,coperture2r=0,coperture12r=0;

    //the common factor of the partition, AND'ing
    for (int i = 0; i < N; i++)
        common_factor[i] = first->binary[i] & second->binary[i];
    common_factor[0] = 1;

    //reducing both the partitions XOR'ing with the common factor
    for (int i = 0; i < N; i++)
        reduced1[i] = common_factor[i] ^ first->binary[i];
    reduced1[0] = 1;
    for (int i = 0; i < N; i++)
        reduced2[i] = common_factor[i] ^ second->binary[i];
    reduced2[0] = 1;

    //the partition product/intersection, OR'ing
    for (int i = 0; i < N; i++)
        product[i] = first->binary[i] | second->binary[i];
    product[0] = 1;

    //intersection of the reduced partitions, OR'ing
    for (int i = 0; i < N; i++)
        product_reduced[i] = reduced1[i] | reduced2[i];
    //  alternatively, one could XOR the unreduced partitions
    //  intersection_reduced[i] =first->binary[i] ^ second->binary[i];
    product_reduced[0] = 1;

    for (int i = 0; i < N; i++) {
        coperture1 += first->binary[i];
        coperture2 += second->binary[i];
        coperture12+=product[i];
        coperture1r+=reduced1[i];
        coperture2r+=reduced2[i];
        coperture12r += product_reduced[i];
    }
    
    res.dist_ham=0;
    for(int i=0;i<N;i++)
        res.dist_ham+= (double) (e1.seq[i] != e2.seq[i]);
    

    //calculating ALL the entropies (3 per non-reduced Rohlin dist, 3 per reduced)
    double
    h1 = entropy_binary_partition(first->binary, N),
            h2 = entropy_binary_partition(second->binary, N),
            hr1 = entropy_binary_partition(reduced1, N),
            hr2 = entropy_binary_partition(reduced2, N),
            h12 = entropy_binary_partition(product, N),
            hr12 = entropy_binary_partition(product_reduced, N);

    //packing the output tuple with the distances
    res.dist_s = 2 * h12 - h1 - h2;
    res.dist_s_r = 2 * hr12 - hr1 - hr2;
    res.dist_top = 2*mylog[coperture12]-mylog[coperture1]-mylog[coperture2] ;
    res.dist_top_r = 2*mylog[coperture12r]-mylog[coperture1r]-mylog[coperture2r];

    //printing the pretty partition graphs if selected
    if (opts.verbose > 2) {
        printf("H(1):\t\t%f\n", h1);
        print_binary_partition(first->binary, N);
        printf("H(2):\t\t%f\n", h2);
        print_binary_partition(second->binary, N);
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

        printf("d(1,2)=\t\t%f\n", res.dist_s);
        printf("dr(1,2)=\t%f\n", res.dist_s_r);
        printf("dt(1,2)=\t%f\n", res.dist_top);
        printf("dtr(1,2)=\t%f\n", res.dist_top_r);
    }

    return (res);
}

dist_tuple general_distance(general_partition *p1, general_partition *p2){
    dist_tuple res;
        
    for(int i=0;i<p1->N;i++)
        product[i]=p1->labels[i]+2*(p2->labels[i]);
    
    general_partition *p12=new general_partition(product,p1->N);
        
    double  h12=p12->entropia_shannon(),
            h1=p1->entropia_shannon(),
            h2=p2->entropia_shannon(),
            t12=p12->entropia_topologica(),
            t1=p1->entropia_topologica(),
            t2=p2->entropia_topologica();
    res.dist_fuzzy=2*h12-h1-h2;
    res.dist_fuzzy_t=2*t12-t1-t2;
    
    //p12->print();
    
    //printf("h12: %f, h1: %f, h2:%f\n",h12,h1,h2);
    //printf("%f\n",2*h12-h1-h2);
    //printf("t12: %f, t1: %f, t2:%f\n",t12,t1,t2);
    //printf("dist: %f\n\n\n",2*t12-t1-t2);
    
    delete p12;
    return(res);
}

