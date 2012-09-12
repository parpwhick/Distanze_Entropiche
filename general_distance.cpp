/**@file general_distance.cpp
 *
 * @author Dawid Crivelli
 *
 * @brief File main(), legge, inizializza, calcola le partizioni e i dati necessari
 * @date 11/01/2012
 *
 *
 * Version: 7.0, 2012/08/04
 * -Complete documentation with Doxygen
 * -Nicer output from time_series() and full command line control
 *
 * Version: 6.3, 2012/07/29
 * -Memory efficience guaranteed
 * -Improvements in findroot and is_equal
 * -Unified entropy calculation
 * -Product of partitions
 * -Cleanup and comments in the code
 *
 * Version: 6.2, 2012/07/13
 * -Improved simulation, with time series output
 * -23x speed improvement over version 6.0!
 * -Using STL classes and methods
 * -Using a more optimized reduction (without common part)
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
#include "strutture.h"
#include "ising_simulation.h"
#include "partizioni.h"
#include "adj_handler.h"
#include "distance.h"

/// Variabile globale con le impostazioni
extern options opts;
/** @brief Tabella lookup dei logaritmi
 *
 * L'utilizzo del lookup aumenta di 6x la velocità del programma!
 * È definita globalmente e calcolata nel main().
 */
double *mylog = 0;

void demo(const adj_struct &topology);

/**
 * @brief Stampa numero medio di atomi nelle partizioni generate
 * @param X Array di partizioni lineari
 */
void print_partition_stats(general_partition *X) {
    label_t min = X[0].N * 20, max = 0;
    double mean = 0, std = 0;
    for (int i = 0; i < opts.n_seq; i++) {
        if (X[i].n < 10)
            printf("Few atoms in partition: %d\n", i);
        min = std::min(X[i].n, min);
        max = std::max(X[i].n, max);
        mean += X[i].n;
        std += (X[i].n)*(X[i].n);

    }
    mean /= opts.n_seq;
    std /= opts.n_seq;
    std = std - mean*mean;

    fprintf(stderr,"Partitions: n. atoms between [%d,%d], ",min,max);
    fprintf(stderr,"avg %.2f with %.2f sites/atom\n",mean,opts.seq_len/mean);

}

/**
 * @brief Funzione main per la manipolazioni di oggetti @ref general_partition. Legge le configurazioni, stampa le statistiche, crea le partizioni, calcola le distanze.
 *
 */
int main(int argc, char** argv) {
    ///Imposta il tipo di partizione
    opts.partition_type = GENERAL_PARTITION;
    ///Interpreta le opzioni da linea di comando
    set_program_options(opts, argc, argv);

    //
    //  ALLOCATION OF MEMORY AND INITIALIZATION
    //
    general_partition *Z = new general_partition[opts.n_seq];

    ///Carica l'opportuna struttura di adiacenza, selezionata da linea di comando
    adj_struct topologia;
    switch (opts.topologia) {
        case(RETICOLO_2D):
            topologia = adiacenza_square_lattice(opts.lato);
            break;
        case(TORO_2D):
            topologia = adiacenza_toroidal_lattice(opts.lato);
            break;
        case(SIERPINSKI):
            topologia = adiacenza_sierpinski(opts.sierpinski_gen);
            break;
        case(FUZZY):
            topologia = adiacenza_fuzzy_line(opts.seq_len);
            break;
        case(LINEARE):
            topologia = adiacenza_simple_line(opts.seq_len);
            break;
        default:
        case(FROM_FILE):
            topologia = adiacenza_from_file(opts.adj_vec_1, opts.adj_vec_2);
            break;
    }
    ///Il volume della partizione è determinato dalla sua struttura di adiacenza
    opts.seq_len = topologia.N;

    ///Stima la quantità di memoria utilizzata, per poter valutare se è sufficiente
    unsigned long memory_estimate =
            +topologia.n_link * 2 * sizeof (int) +topologia.N * sizeof (int) //struct adiacenza
            +opts.seq_len * sizeof (int) //caricamento sequenza
            +opts.seq_len * sizeof (double) //logaritmo
            +opts.seq_len * sizeof (product_t) * opts.threads; //temp product
    if(opts.letto_da != SIMULATION)
        memory_estimate+=
             opts.seq_len * opts.n_seq * sizeof (label_t) * ( 2 +1) //partizioni: fisso + max_variabile/3
            +opts.n_seq * opts.n_seq * __builtin_popcount(opts.da_calcolare) * sizeof(double); //per l'output
    else
        memory_estimate+=
             opts.seq_len * 2 * sizeof (label_t) * ( 2 +1) //partizioni: fisso + max_variabile/3
            + sizeof(ising_simulation::config_t) * (topologia.N + topologia.n_link) ;
    memory_estimate >>= 20;
    fprintf(stderr,"Estimated memory usage: %lu MB\n",memory_estimate+1);

    //logarithm lookup table, 6x program speedup
    int lunghezza = std::max(opts.seq_len, opts.n_seq) + 10;
    mylog = new double[ lunghezza ];
    for (int i = 1; i < lunghezza; i++)
        mylog[i] = log(i);
    mylog[0] = 0;

    //
    // DEMO
    //
    if (opts.demo) {
        demo(topologia);
        return 0;
    }
    //
    // SIMULATION case
    //
    /**In caso di configurazioni generate tramite simulazione, non calcola la matrice completa delle distanze,
     * ma la time_series() tra configurazioni e partizioni successive nel tempo.
     */
    if(opts.letto_da == SIMULATION){
        time_series(topologia);
        return(0);
    }

    //
    // RANDOM or FROM_FILE
    //
    int *num_buffer = new int[opts.seq_len];
    for (int i = 0; i < opts.n_seq; i++) {

        if (opts.letto_da == FROM_FILE) {
            ///Se le configurazioni sono da leggere da file, ne carica una alla volta e ne calcola la partizione
            int result = load_config(opts, num_buffer);
            if (!result)
                break;
        }
        if (opts.letto_da == RANDOM)
            ///In caso di configurazioni random, ne crea una per volta
            generate_next_sequence(num_buffer);

        if (opts.verbose)
            fprintf(stderr, "Loaded sequence %d, analysing\n", i + 1);
        ///Generazione della partizione
        Z[i].from_configuration(num_buffer, topologia);
    }

    printf("Loaded %d sequences long %d\n", opts.n_seq, opts.seq_len);
    print_partition_stats(Z);
    printf("\n");

    //
    //  DISTANCE MEASUREMENTS
    //
    if (opts.distance == false)
        return(0);

    calcola_matrice_distanze(Z);
    //
    //  PROGRAM EXIT
    //
    delete []Z;
    delete []mylog;
    delete []num_buffer;
    return 0;
}

#include <iostream>
#include <iomanip>
#include <ctime>
using namespace std;
void demo(const adj_struct &topology) {
    general_partition A;
    general_partition B;
    general_partition comune, prodotto;

    opts.seq_len = topology.N;
    vector<int> num_buffer(opts.seq_len);

    while (1) {
        generate_next_sequence(num_buffer.data());
        A.from_configuration(num_buffer.data(), topology);
        generate_next_sequence(num_buffer.data());
        B.from_configuration(num_buffer.data(), topology);

        comune.common_subfactor(B, A);
        prodotto.product(A, B);
        //if(comune.n > 1)
        if (!(A == B))
            break;
    }

    A.consistency_check();
    B.consistency_check();
    prodotto.consistency_check();
    comune.consistency_check();

    cout.setf(ios::boolalpha);
    cout << endl
            << "Demo, using partitions long " << opts.seq_len << endl
            << "A,B - partitions" << endl
            << "A^B - common subfactor" << endl
            << "AvB - product" << endl
            << endl
            << "A   == B  : " << (A == B) << "\n"
            << "A   <= B  : " << (A <= B) << "\n"
            << "A^B <= A  : " << (comune <= A) << "\n"
            << "A^B <= B  : " << (comune <= B) << "\n"
            << "A   <= AvB: " << (A <= prodotto) << "\n"
            << "B   <= AvB: " << (B <= prodotto) << "\n"
            << endl;

    general_partition temp1;
    temp1.common_subfactor(prodotto, B);
    cout << "(AvB)^B == B: " << (temp1 == B) << "\n";
    temp1.product(comune, A);
    cout << "(A^B)vA == A: " << (temp1 == A) << "\n";
    cout << endl;

    std::clock_t start;
    double time_diff = 0;
    int RETRIES = 0;

    start = std::clock();
    while (time_diff < 250) {
        temp1.product(A, B);
        RETRIES++;
        time_diff = (std::clock() - start) * 1000. / (double) CLOCKS_PER_SEC;
    }
    RETRIES *= 4;
    RETRIES += 2;

    start = std::clock();
    for (int k = 0; k < RETRIES; k++)
        temp1.product(A, B);
    time_diff = (std::clock() - start) * 1000. / (double) CLOCKS_PER_SEC / RETRIES;
    fprintf(stdout, "Product      in %.3g ms\n", time_diff);

    /*start = std::clock();
    for(int k=0; k< RETRIES; k++)
        temp1.quick_product(A,B);
    time_diff = (std::clock() - start) / (double) CLOCKS_PER_SEC / RETRIES;
    fprintf(stdout, "Quick prod   in %.3e seconds\n", time_diff);*/

    start = std::clock();
    for (int k = 0; k < RETRIES; k++)
        temp1.common_subfactor(A, B);
    time_diff = (std::clock() - start) * 1000. / (double) CLOCKS_PER_SEC / RETRIES;
    fprintf(stdout, "Intersection in %.3g ms\n", time_diff);

    start = std::clock();
    for (int k = 0; k < RETRIES; k++)
        temp1.reduce(A, B);
    time_diff = (std::clock() - start) * 1000. / (double) CLOCKS_PER_SEC / RETRIES;
    fprintf(stdout, "Reduction    in %.3g ms\n", time_diff);

    ::distance dist(A.N);
    opts.da_calcolare = SHAN;
    start = std::clock();
    for (int k = 0; k < RETRIES; k++)
        dist(A, B);
    time_diff = (std::clock() - start) * 1000. / (double) CLOCKS_PER_SEC / RETRIES;
    fprintf(stdout, "Dist(short)  in %.3g ms\n", time_diff);

    opts.da_calcolare = SHAN | RID;
    start = std::clock();
    for (int k = 0; k < RETRIES; k++)
        dist(A, B);
    time_diff = (std::clock() - start) * 1000. / (double) CLOCKS_PER_SEC / RETRIES;
    fprintf(stdout, "Dist(+red)   in %.3g ms\n", time_diff);
}
