
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>

using namespace std;
int verbose = 0;
double *mylog;

struct dist_tuple {
    double dist;
    double dist_r;
    double dist_top;
    double dist_top_r;
};

int *match;
int *reduced1;
int *reduced2;
int *intersection;
int *intersection_reduced;

class partition {
public:
    //array of atom starting positions {1,4,5,...}
    int *atom_positions;
    //binary atom beginning codes {1,0,0,0,1,1,...}
    int *binary;
    //number of atoms found
    int n;
    //total length of the partition
    int N;

    double entropia_shannon();
    int entropia_topologica() {
        return n;
    }
    void print();
    partition() {
        atom_positions = 0;
        binary = 0;
        n = 0;
        N = 0;
    }

};


class entry {
public:
    //identifier of the sequence
    char* name;
    //string of symbols making the sequence
    char* seq;
    //length of the sequence
    int n;
    //partitioned sequence
    partition Z;

    void make_partition();
    void print();
};

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
        match[i] = first.binary[i] & second.binary[i];
    match[0] = 1;

    //reducing both the partitions XOR'ing with the common factor
    for (int i = 0; i < N; i++)
        reduced1[i] = match[i] ^ first.binary[i];
    reduced1[0] = 1;
    for (int i = 0; i < N; i++)
        reduced2[i] = match[i] ^ second.binary[i];
    reduced2[0] = 1;

    //the partition product/intersection, OR'ing
    for (int i = 0; i < N; i++)
        intersection[i] = first.binary[i] | second.binary[i];
    intersection[0] = 1;

    //intersection of the reduced partitions, OR'ing
    for (int i = 0; i < N; i++)
        intersection_reduced[i] = reduced1[i] | reduced2[i];
    //  alternatively, one could XOR the unreduced partitions
    //  intersection_reduced[i] =first.binary[i] ^ second.binary[i];
    intersection_reduced[0] = 1;

    for (int i = 0; i < N; i++) {
        coperture1 += first.binary[i];
        coperture2 += second.binary[i];
        coperture12+=intersection[i];
        coperture1r+=reduced1[i];
        coperture2r+=reduced2[i];
        coperture12r += intersection_reduced[i];
    }

    //calculating ALL the entropies (3 per non-reduced Rohlin dist, 3 per reduced)
    double
    h1 = entropy_binary_partition(first.binary, N),
            h2 = entropy_binary_partition(second.binary, N),
            hr1 = entropy_binary_partition(reduced1, N),
            hr2 = entropy_binary_partition(reduced2, N),
            h12 = entropy_binary_partition(intersection, N),
            hr12 = entropy_binary_partition(intersection_reduced, N);

    //packing the output tuple with the distances
    res.dist = 2 * h12 - h1 - h2;
    res.dist_r = 2 * hr12 - hr1 - hr2;
    res.dist_top = 2*mylog[coperture12]-mylog[coperture1]-mylog[coperture2] ;
    res.dist_top_r = 2*mylog[coperture12r]-mylog[coperture1r]-mylog[coperture2r];

    //printing the pretty partition graphs if selected
    if (verbose > 2) {
        printf("H(1):\t\t%f\n", h1);
        print_binary_partition(first.binary, N);
        printf("H(2):\t\t%f\n", h2);
        print_binary_partition(second.binary, N);
        if (verbose > 3) {
            printf("Hr(1):\t\t%f\n", hr1);
            print_binary_partition(reduced1, N);
            printf("Hr(2):\t\t%f\n", hr2);
            print_binary_partition(reduced2, N);
        }
        printf("H(12):\t\t%f\n", h12);
        print_binary_partition(intersection, N);
        printf("Hr(12):\t\t%f\n", hr12);
        print_binary_partition(intersection_reduced, N);
        printf("\n");

        printf("d(1,2)=\t\t%f\n", res.dist);
        printf("dr(1,2)=\t%f\n", res.dist_r);
        printf("dt(1,2)=\t%f\n", res.dist_top);
        printf("dtr(1,2)=\t%f\n\n\n", res.dist_top_r);
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
    fprintf(stderr, "%s: %s\n", name, seq);
    //printing of the whole partition
    Z.print();
    //adding the entropy information
    fprintf(stderr, "Shannon %f, Topological %d\n\n", Z.entropia_shannon(),
            Z.entropia_topologica());
}

/*
 * 
 */
int main(int argc, char** argv) {

    int i, j;
    int seq_len = 100;
    int n_seq = 1000;
    double *distances;
    double *distances_r;
    double *dist_top;
    double *dist_top_r;
    FILE *out;
    dist_tuple res;

    //
    //  ALLOCATION OF MEMORY AN INITIALIZATION
    //
    entry *entries = new entry[n_seq];

    //arrays needed to intersect the partitions
    match = new int[seq_len];
    reduced1 = new int[seq_len];
    reduced2 = new int[seq_len];
    intersection = new int[seq_len];
    intersection_reduced = new int[seq_len];

    //logarithm lookup table, 6x program speedup
    mylog = new double[seq_len + 10];
    for (i = 0; i < seq_len + 10; i++)
        mylog[i] = log(i);

    //random number init
    srand(37337);

    //distance matrix initialization and zeroing
    distances = new double[n_seq * n_seq];
    distances_r = new double[n_seq * n_seq];
    dist_top = new double[n_seq * n_seq];
    dist_top_r = new double[n_seq * n_seq];
    for (i = 0; i < n_seq * n_seq; i++) {
        distances[i] = 0;
        distances_r[i] = 0;
        dist_top[i]=0;
        dist_top_r[i]=0;
    }

    
    
    
    //
    //  CREATION OF RANDOM SEQUENCES (as entries)
    //
    for (i = 0; i < n_seq; i++) {
        char *name = new char[10];
        snprintf(name, 10, "seq%03d", i + 1);
        entries[i].name = name;
        char *seq = new char[seq_len + 1];
        for (j = 0; j < seq_len; j++)
            seq[j] = rand() % 2 + 'A';
        seq[seq_len + 1] = 0;
        entries[i].seq = seq;
        entries[i].n = seq_len;
    }

    //
    //  ANALISYS OF LOADED SEQUENCES
    //
    for (i = 0; i < n_seq; i++) {
        entries[i].make_partition();
        if (verbose > 0)
            entries[i].print();
    }

    //
    //  DISTANCE MEASUREMENTS
    //
    for (i = 0; i < n_seq; i++)
        for (j = i + 1; j < n_seq; j++) {
            res = intersect(entries[i].Z, entries[j].Z);
            distances[i * n_seq + j] = res.dist;
            distances_r[i * n_seq + j] = res.dist_r;
            dist_top[i*n_seq+j]=res.dist_top;
            dist_top_r[i*n_seq+j]=res.dist_top_r;
        }


    //
    //  WRITING THE OUTPUT
    //
    out = fopen("output-dist.bin", "w");
    i = fwrite(distances, sizeof (double), n_seq*n_seq, out);
    printf("Written %d records\n", i);
    fclose(out);

    out = fopen("output-distr.bin", "w");
    i = fwrite(distances_r, sizeof (double), n_seq*n_seq, out);
    printf("Written %d records\n", i);
    fclose(out);

    out = fopen("output-top.bin", "w");
    i = fwrite(dist_top, sizeof (double), n_seq*n_seq, out);
    printf("Written %d records\n", i);
    fclose(out);
    
    out = fopen("output-topr.bin", "w");
    i = fwrite(dist_top_r, sizeof (double), n_seq*n_seq, out);
    printf("Written %d records\n", i);
    fclose(out);
    //
    //  PROGRAM EXIT
    //
    return 0;
}


