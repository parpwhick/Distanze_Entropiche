/** \file distance.cpp
 * @brief Definizione dei metodi per la classe @ref distance e l'essenziale @ref calcola_matrice_distanze().
 * Definisce anche @ref entropy_binary_partition()
 */
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <ctime>
#include <vector>
#include <algorithm>
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "strutture.h"
#include "partizioni.h"
#include "distance.h"

extern options opts;
extern double *mylog;

/* Prealloca i vettori temporanei utilizzati per il calcolo delle entropie
 */
void distance::allocate(int n) {
    if (opts.partition_type == LINEAR_PARTITION) {
        common_factor.resize(n);
        reduced1.resize(n);
        reduced2.resize(n);
        product_reduced.resize(n);
        binary_product.resize(n);
    } else {
        label_index.resize(n);
    }
}

distance::distance(int n) {
    dist_ham=0;
    dist_shan=0;
    dist_shan_r=0;
    dist_top=0;
    dist_top_r=0;
    N = n;
    allocate(N);
}

distance::distance(const distance &d1) {
    N = d1.N;
    allocate(N);
}

/** Calcola tutte le distanze tra le due partizioni. Per fare questo prima calcola
 *  le partizioni ridotte, successivamente ne esegue la distanza. Infine calcola
 *  la distanza tra le partizioni non ridotte.
 */
void distance::dist(const general_partition& e1, const general_partition& e2) {
    if (opts.da_calcolare & RID) {
        if (opts.riduzione == COMUNE) {
            partizione_comune.common_subfactor(e1, e2);
            ridotto1.reduce(e1, partizione_comune);
            ridotto2.reduce(e2, partizione_comune);
        } else if (opts.riduzione == DIRETTA) {
            ridotto1.reduce(e1, e2);
            ridotto2.reduce(e2, e1);
        }
        //stampa delle etichette in modo testuale, eventualmente da rimuovere
        if (opts.verbose > 2) {
            label_t quanto = std::min(e1.N, (label_t) 50);
            print_array(&partizione_comune.labels[0], quanto, "lbls comune");
            print_array(&e1.labels[0], quanto, "lbls e1    ");
            print_array(&ridotto1.labels[0], quanto, "lbls ridot1");
            print_array(&e2.labels[0], quanto, "lbls e2    ");
            print_array(&ridotto2.labels[0], quanto, "lbls ridot2");
        }

        calc_distance(ridotto1, ridotto2);
        dist_shan_r = dist_shan;
        dist_top_r = dist_top;
        //stampa grafici delle partizioni ridotte e nonridotte, forse da rimuovere
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
        calc_distance(e1, e2);
}

/**@brief Stampa rappresentazione di una partizione binaria: |...|...||...|...
 *
 * @param p Array con la partizione
 * @param N Lunghezza dell'array
 */
void print_binary_partition(int*p, int N) {
    for (int i = 0; i < N; i++)
        printf("%c", (p[i]) ? '|' : '.');
    printf("\n");
}

/** @brief Calcolo dell'entropia di un vettore di 1 e 0, in cui 1 indica un nuovo atomo
 *
 * La funzione molto semplice per il calcolo dell'entropia di una partizione
 * lineare rappresentata in modo binario, in cui l'inizio di un nuovo atomo e' indicato da "1", a differenza
 * delle simili funzioni utilizzate altrove
 * @param p Vettore con la partizione binaria
 * @param N Lunghezza del vettore
 * @return H Entropia di Shannon
 */
template <typename T>
inline double entropy_binary_partition(const std::vector<T> &p, int N) {
    int label_count = 0;
    double H = 0;
    int mu;
    int begin;

    //the first position always starts an atom
    begin = 0;
    label_count = 1;

    for (label_t i = 1; i < N; i++) {
        //whenever we find a new atom
        if (p[i]) {
            //the closed (old)atom's length is calculated
            mu = i - begin;
            //the new one is ready to go
            label_count++;
            begin = i;
            //we add the entropy, with the trick mu>0 and when mu=1 the log is 0
            if (mu > 1)
                H += (double) mu * mylog[mu];
        }
    }
    //the last one, so it's not left hanging
    mu = N - begin;
    H += mu * mylog[mu];

    //normalize the result
    H = -H / N + mylog[N];

    return H;
}

template <typename T>
void distance::hamming_distance(const T* seq1, const T* seq2) {
    this->dist_ham = 0;
    for (int i = 0; i < N; i++)
        this->dist_ham += (double) (seq1[i] != seq2[i]);
}
template void distance::hamming_distance(const char*, const char*);

/**
 * Calcola le 4 possibili distanze tra partizioni lineari - Shannon, Shannon ridotto,
 * topologica, topologica ridotta.
 *
 * I metodi sono estremamente efficienti, corrispondenti a operazioni bitwise e
 * calcolo di entropie (che corrispondono alla maggior parte del tempo utilizzato).
 *
 * @param first Prima partizione lineare
 * @param second Seconda partizione lineare
 * @return Le distanze sono scritte nelle variabili interne dell'oggetto @c distance
 */
void distance::dist(const linear_partition &first, const linear_partition &second) {

    int N = first.N;
    bool ridotta = opts.da_calcolare & RID;
    //the partition product/intersection, OR'ing
    for (int i = 0; i < N; i++)
        binary_product[i] = first.binary[i] | second.binary[i];
    binary_product[0] = 1;

    //calculating ALL the entropies (3 per non-reduced Rohlin dist, 3 per reduced)
    double
    h1 = first.entropia_shannon,
            h2 = second.entropia_shannon,
            h12 = entropy_binary_partition(binary_product, N);

    this->dist_shan = 2 * h12 - h1 - h2;


    //DISTANZA TOPOLOGICA
    //per "coperture" si intende il numero di atomi di una partizione
    int coperture1 = 0, coperture2 = 0, coperture12 = 0;
    for (int i = 0; i < N; i++) {
        coperture1 += first.binary[i];
        coperture2 += second.binary[i];
        coperture12 += binary_product[i];
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
    //  product_reduced[i] =first.binary[i] ^ second.binary[i];
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

    //printf("comune semplice: %d, r1: %d/%d, r2: %d/%d, prod: %d\n",
    //coperture_common, coperture1r, coperture1, coperture2r, coperture2, coperture12r);
}

/**
 * Calcola la distanza tra due partizioni, utilizzando solamente i labels. Il metodo
 * utilizzato consiste nel calcolare i label della partizione prodotto ed analizzare
 * la loro entropia con il metodo distruttivo del sort.
 * L'algoritmo e' incredibilmente
 * piu' veloce di qualunque metodo alternativo -- sort accede alla memoria di un vettore
 * in maniera ottimale, in modo impossibile da battere usando mappe ordinate di coppie
 * o simili trucchi.
 * @param p1 Prima fattore
 * @param p2 Secondo fattore
 */
void distance::calc_distance(const general_partition &p1, const general_partition &p2) {
    label_t count = 0;
    label_t old_count = 0;
    int fiddle = 0;
    for (label_t atom_index = 0; atom_index < p1.n; atom_index++) {
        fiddle = (fiddle) ? 0 : (1<<31);
        for (Iter_t ii = p1.begin(atom_index); ii != p1.end(); ii++)
            label_index[count++] = fiddle | p2.labels[*ii];

        std::sort(&label_index[old_count], &label_index[count]);
        old_count = count;
    }

    std::pair<double, int> entropie = ordered_vector_entropy(label_index.data(), N);
    int n = entropie.second;
    double H = entropie.first;
    
    double h1 = p1.entropia_shannon,
            h2 = p2.entropia_shannon,
            t1 = p1.entropia_topologica,
            t2 = p2.entropia_topologica;
    this->dist_shan = 2 * H - h1 - h2;
    this->dist_top = 2 * mylog[n] - t1 - t2;

    if (opts.graphics && (opts.topologia == RETICOLO_2D)) {
        static int imagecount = 0;
        char filename[255];
        imagecount++;
        sprintf(filename, "prodotto%03d.ppm", imagecount);
        ppmout2(&p1.labels[0], &p2.labels[0], opts.lato, filename);
        imagecount++;
        sprintf(filename, "prodotto%03d.ppm", imagecount);
        ppmout(&product[0], opts.lato, filename);
    }
}

/**@brief
 * Semplice funzione template per scrivere brevemente in molti file
 * @param where Nome del file in cui scrivere
 * @param what Vettore contenente i dati da scrivere
 * @return
 */
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

/** @brief Dato un vettore di partizioni, ne scrive la matrice completa delle distanze.
 *
 * Dato un vettore di partizioni, semplici o generali, la funzione alloca la memoria necessaria,
 * gli oggetti @ref distance e calcola in modo multithreaded la matrice delle distanze tra le
 * partizioni considerate. Un'altra utile caratteristica e' la stima del tempo richiesto in
 * real-time. Alla fine scrive (se necessario) i valori calcolati.
 *
 * @param X Vettore di partizioni (di qualunque tipo)
 */
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
    if (da_calcolare & SHAN) dist_shan.resize(dim);
    if (da_calcolare & RID) dist_shan_r.resize(dim);
    if (da_calcolare & TOP) dist_top.resize(dim);
    if (da_calcolare & RID_TOP) dist_top_r.resize(dim);
    if (da_calcolare & HAMM) dist_ham.resize(dim);

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
                d.dist(X[i], X[j]);
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
