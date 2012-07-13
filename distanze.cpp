#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <ctime>
#include <vector>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "strutture.h"

extern options opts;
extern double *mylog;

void distance::allocate(int n) {
    if (opts.partition_type == LINEAR_PARTITION) {
        common_factor.reserve(n);
        reduced1.reserve(n);
        reduced2.reserve(n);
        product_reduced.reserve(n);
    }
    product.reserve(n);
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
    N = 0;
}

void distance::fill(const linear_partition& e1, const linear_partition& e2) {
    this->dist(e1, e2);
}

void distance::fill(const general_partition& e1, const general_partition& e2) {
    

    if (opts.da_calcolare & RID) {
        //partizione_comune.linear_intersection(e1, e2);
        //partizione_comune.trivial(e1.N);
        //ridotto1.reduce(e1, partizione_comune);
        //ridotto2.reduce(e2, partizione_comune);


        if (opts.verbose > 2) {
            label_t quanto = std::min(e1.N, (label_t) 50);
            print_array(&partizione_comune.labels[0], quanto, "lbls comune");
            print_array(&e1.labels[0], quanto, "lbls e1    ");
            print_array(&ridotto1.labels[0], quanto, "lbls ridot1");
            print_array(&e2.labels[0], quanto, "lbls e2    ");
            print_array(&ridotto2.labels[0], quanto, "lbls ridot2");

        }
                ridotto1.reduce(e1, e2);
                ridotto2.reduce(e2, e1);

        dist(ridotto1, ridotto2);
        dist_shan_r = dist_shan;
        dist_top_r = dist_top;
        //printf("comune generale: %d, r1: %d/%d, r2: %d/%d, prod: %d\n",partizione_comune.n, ridotto1.n, e1.n, ridotto2.n,e2.n,(int)dist_fuzzy_t);
        if (opts.graphics && (opts.topologia == RETICOLO_2D)) {
            static int imagecount = 0;
            char filename[255];
            imagecount++;
            sprintf(filename, "ridotto%03d.ppm", imagecount);
            ppmout2(&e1.labels[0], &e2.labels[0], opts.lato, filename);
            imagecount++;
            sprintf(filename, "ridotto%03d.ppm", imagecount);
            ppmout2(&ridotto1.labels[0], &ridotto2.labels[0], opts.lato, filename);
        }
    }

    if (opts.da_calcolare & SHAN)
        dist(e1, e2);
}

void print_binary_partition(int*p, int N) {
    // |...|...||...|...
    for (int i = 0; i < N; i++)
        printf("%c", (p[i]) ? '|' : '.');
    printf("\n");
}

template <typename T>
double entropy_binary_partition(const T *p, int N) {
    //p: binary partition
    //n: total length
    int i = 0;
    int mu;
    int begin;
    double H = 0;

    //the first position always starts an atom
    begin = 0;
    for (i = 1; i < N; i++) {
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
    mu = N - begin;
    if (mu > 1) H += mu * mylog[mu];

    //proper entropy normalization
    H = -H / N + mylog[N];
    return (H);
}
template double entropy_binary_partition(const int *, int);

template <typename T>
inline double entropy_binary_partition(const std::vector<T> &p, int N) {
    return entropy_binary_partition(&p[0],N);    
}

template <typename T>
void distance::hamming_distance(const T* seq1, const T* seq2) {
    this->dist_ham = 0;
    for (int i = 0; i < N; i++)
        this->dist_ham += (double) (seq1[i] != seq2[i]);
}
template void distance::hamming_distance(const char*, const char*);



//function giving the distance between 2 partitions, by intersecting 
//and calculating the relevant entropy

void distance::dist(const linear_partition &first, const linear_partition &second) {

    int N = first.N;
    bool ridotta = opts.da_calcolare & RID;
    //the partition product/intersection, OR'ing
    for (int i = 0; i < N; i++)
        product[i] = first.binary[i] | second.binary[i];
    product[0] = 1;

    //calculating ALL the entropies (3 per non-reduced Rohlin dist, 3 per reduced)
    double
    h1 = first.entropia_shannon,
            h2 = second.entropia_shannon,
            h12 = entropy_binary_partition(product, N);

    this->dist_shan = 2 * h12 - h1 - h2;


    //DISTANZA TOPOLOGICA
    int coperture1 = 0, coperture2 = 0, coperture12 = 0;
    for (int i = 0; i < N; i++) {
        coperture1 += first.binary[i];
        coperture2 += second.binary[i];
        coperture12 += product[i];
    }
    this->dist_top = 2 * mylog[coperture12] - mylog[coperture1] - mylog[coperture2];


    if (!ridotta)
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

    double hr1 = entropy_binary_partition(reduced1, N),
            hr2 = entropy_binary_partition(reduced2, N),
            hr12 = entropy_binary_partition(product_reduced, N);
    this->dist_shan_r = 2 * hr12 - hr1 - hr2;


    //DISTANZA TOPOLOGICA RIDOTTTA
    int coperture1r = 0, coperture2r = 0, coperture12r = 0, coperture_common = 0;
    for (int i = 0; i < N; i++) {
        coperture1r += reduced1[i];
        coperture2r += reduced2[i];
        coperture12r += product_reduced[i];
        coperture_common += common_factor[i];
    }
    this->dist_top_r = 2 * mylog[coperture12r] - mylog[coperture1r] - mylog[coperture2r];

    //printf("comune semplice: %d, r1: %d/%d, r2: %d/%d, prod: %d\n",coperture_common, coperture1r, coperture1, coperture2r, coperture2, coperture12r);
}

void distance::dist(const general_partition &p1, const general_partition &p2) {
    int i;
    int label_count = 0;
    double H = 0;
    int mu;
    int begin;

    for (i = 0; i < N; i++) {
        uint64_t temp1 = p1.labels[i];
        uint64_t temp2 = p2.labels[i];

        product[i] = (temp1 << 32) | temp2;
    }

    static int imagecount = 0;
    char filename[255];
    if (opts.graphics && (opts.topologia == RETICOLO_2D)) {
        imagecount++;
        sprintf(filename, "prodotto%03d.ppm", imagecount);
        ppmout2(&p1.labels[0], &p2.labels[0], opts.lato, filename);
        imagecount++;
        sprintf(filename, "prodotto%03d.ppm", imagecount);
        ppmout(&product[0], opts.lato, filename);
    }

    std::sort(product.begin(), product.end());

    //the first position always starts an atom
    begin = 0;
    label_count = 1;
    uint64_t old_val = product[0];
    for (i = 1; i < N; i++) {
        //whenever we find a different value (a new atom)
        if (product[i] != old_val) {
            //a new atom starts
            label_count++;
            //the closed (old)atom's length is calculated
            mu = i - begin;
            //the new one is ready to go
            begin = i;
            //cache the new label to check
            old_val = product[i];
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
    this->dist_shan = 2 * H - h1 - h2;
    this->dist_top = 2 * mylog[label_count] - t1 - t2;


}

int WRITE(const char *where, const std::vector<double> & what) {
    FILE *out = fopen(where, "wb");
    int expected = opts.n_seq * opts.n_seq;
    int bytes_written = fwrite(&what[0], sizeof (double), expected, out);
    fclose(out);
    if (bytes_written == expected)
        return (1);
    else {
        printf("Expected to write %d bytes, %d instead\n", expected, bytes_written);
        return (0);
    }
}

template <typename partition_t>
void calcola_matrice_distanze(const partition_t* X) {

    int &da_calcolare = opts.da_calcolare;
    const int dim=opts.n_seq * opts.n_seq;

    std::vector<double> dist_shan;
    std::vector<double> dist_shan_r;
    std::vector<double> dist_top;
    std::vector<double> dist_top_r;
    std::vector<double> dist_ham;
    
    //distance matrix allocation and zeroing,if needed
    if (da_calcolare & SHAN) dist_shan.reserve(dim);
    if (da_calcolare & RID) dist_shan_r.reserve(dim);
    if (da_calcolare & TOP) dist_top.reserve(dim);
    if (da_calcolare & RID_TOP) dist_top_r.reserve(dim);
    if (da_calcolare & HAMM) dist_ham.reserve(dim);

    fprintf(stderr, "Calculating distance matrix\n");

    distance d(opts.seq_len);

    std::clock_t start = std::clock();
    double time_diff;
    double completed_ratio;

#pragma omp parallel for firstprivate(d) schedule(dynamic) num_threads(opts.threads)
    for (int i = 0; i < opts.n_seq; i++) {
        for (int j = i + 1; j < opts.n_seq; j++) {
            int index = i * opts.n_seq + j;

            if (da_calcolare & (SHAN | RID)) {
                d.fill(X[i], X[j]);
                if (da_calcolare & SHAN) dist_shan[index] = d.dist_shan;
                if (da_calcolare & RID) dist_shan_r[index] = d.dist_shan_r;
                if (da_calcolare & TOP) dist_top[index] = d.dist_top;
                if (da_calcolare & RID_TOP) dist_top_r[index] = d.dist_top_r;
            }

            //            if (da_calcolare & HAMM) {
            //                d.hamming_distance(char_entries[i].c_str(), char_entries[j].c_str());
            //                dist_ham[index] = d.dist_ham;
            //            }

        }

#ifdef _OPENMP
        int this_thread = omp_get_thread_num();
        if (this_thread)
            continue;
        double time_ratio = omp_get_num_threads();
#else
        double time_ratio = 1.0;
#endif
        fprintf(stderr, "\r");
        time_diff = (std::clock() - start) / (double) CLOCKS_PER_SEC / time_ratio;
        int k = i + 1;
        completed_ratio = (2 * k * opts.n_seq - k * k - k + 0.0) / (opts.n_seq * (opts.n_seq - 1));
        fprintf(stderr, "%.1f%% done, ETA %.0fs    ",
                completed_ratio * 100, ceil(time_diff * (1 / completed_ratio - 1)));
        fflush(stderr);

    }
    time_diff = (std::clock() - start) / (double) CLOCKS_PER_SEC;
    fprintf(stderr, "\r100%% done in %.1f seconds of CPU time\n", time_diff);

    //
    //  WRITING THE OUTPUT
    //
    if (opts.write == true) {
        int count = 0;
        if (da_calcolare & SHAN)
            count += WRITE("output-shan.bin", dist_shan);

        if (da_calcolare & RID)
            count += WRITE("output-shan_r.bin", dist_shan_r);

        if (da_calcolare & TOP)
            count += WRITE("output-top.bin", dist_top);

        if (da_calcolare & RID_TOP)
            count += WRITE("output-top_r.bin", dist_top_r);
        
        if (da_calcolare & HAMM)
            count += WRITE("output-hamm.bin", dist_ham);

        //if (opts.verbose)
        fprintf(stderr, "Written %dx distance matrix\n", count);
    }

}
template void calcola_matrice_distanze(const linear_partition *X);
template void calcola_matrice_distanze(const general_partition *X);