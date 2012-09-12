/**
 * @file   strutture.h
 * @author Dawid Crivelli
 * @date January 13, 2012, 10:51 PM
 *
 * @brief Definisce le opzioni globali (@ref options) e gli enum che caratterizzano le opzioni
 */

#ifndef STRUTTURE_H
#define	STRUTTURE_H

#include <vector>
#include <string>
#include <stdint.h>

///Le etichette delle partizioni sono interi a 32bit con segno (per avere anche etichetta -1). Questo limita a 2^32 il volume delle partizioni.
typedef int32_t label_t;
///Il prodotto è una coppia di labels, quindi è rappresentabile come un intero senza segno a 64bit. Nel caso di labels a 64bit, bisogna usare coppie std::pair<label_t,label_t>
typedef uint64_t product_t;

///Tipo di update per la simulazione
enum simulation_t {
    ///Metropolis
    METROPOLIS,
    ///Microcanonico eventualmente con bordi termostatati
    MICROCANONICAL
};

///Definisce i tipi di partizione usati, il che comporta una diversa scelta di algoritmi
enum part_type {
    ///partizioni generali
    GENERAL_PARTITION,
    ///partizioni lineari rappresentate in modo binario
    LINEAR_PARTITION
};

///Definisce i tipi di distanze da calcolare, usato in options::da_calcolare
enum DISTANZE {
    ///Shannon
    SHAN = 1,
    ///topologica
    TOP = 1 << 1,
    ///Shannon ridotta
    RID = 1 << 2,
    ///topologica ridotta
    RID_TOP = 1 << 3,
    ///Hamming
    HAMM = 1 << 20
};

///Tipi di sorgente per le configurazioni e i tipi di topologie supportate
enum source {
    ///sequenze lineari senza salto
    LINEARE,
    ///reticolo quadrato con C.C periodiche toroidali
    TORO_2D,
    ///reticolo quadrato con bordi aperti
    RETICOLO_2D,
    ///triangolo di Sierpinski di generazione options::sierpinski_gen
    SIERPINSKI,
    ///sequenze lineari con salto (ovvero molti primi vicini)
    FUZZY,
    ///configurazioni generate a random
    RANDOM,
    ///configurazioni o adiacenza letta da file
    FROM_FILE,
    ///configurazioni generate via evoluzione temporale
    SIMULATION,
};

///Tipo di riduzione da usare
enum red_strategy{
    ///tramite partizione comune
    COMUNE,
    ///eliminando atomi simili tra le due partizioni
    DIRETTA
};

///Struttura globale di tutte le opzioni del programma
typedef class {
public:
    ///Lunghezza delle partizioni
    int seq_len;
    ///Numero di configurazioni/partizioni richieste
    int n_seq;
    ///Lunghezza del lato, nel caso di un reticolo quadrato bidimensionale
    int lato;
    ///Generazione del triangono di Sierpinski richiesta per l'adiacenza
    int sierpinski_gen;
    ///Cardinalita' dell'alfabeto utilizzato per la generazione di configurazioni random [default: 2]
    int n_symbols;
    ///Tipo di partizione: lineare o generale
    part_type partition_type;
    ///Parametro epsilon per la riduzione [default: 0]
    int epsilon;

    ///Informazioni sulla generazione delle configurazioni
    source letto_da;
    ///Topologia delle configurazioni richieste
    source topologia;
    ///Nome del file da cui leggere le configurazioni
    char state_filename[255];
    ///Nome del file con le righe della matrice di adiacenza
    char adj_vec_1[255];
    ///Nome del file con le colonne della matrice di adiacenza
    char adj_vec_2[255];
    ///Per la topologia della retta con salto, indica il salto massimo
    int fuzzy;
    ///Tipo di riduzione da usare: partizione comune o diretta [default]
    red_strategy riduzione;

    ///Bitmap delle distanze da calcolare [default: tutte]
    int da_calcolare;
    ///Scrivere i risultati in file? [default: true]
    bool write;
    ///Calcolare le distanze? [default: true]
    bool distance;
    ///Numero dei threads richiesti [default: numero di processori]
    int threads;
    ///Stampare disegni in .pbm per tutte le operazioni svolte
    bool graphics;

    ///Seed del generatore di numeri casuali [default: casuale]
    int seed;
    ///Dinamica Metropolis o microcanonica [default: microcanonica]
    simulation_t simulation_type;
    ///Numero di sweeps in un intervallo temporale [default: 1]
    int sweeps;
    ///Numero di istanti temporali da saltare inizialmente
    int skip;
    ///Temperatura inversa per metropolis [default: 0.45]
    double beta;
    ///Massima energia per link nel caso di distribuzione random uniforme
    int max_link_energy;

    ///Grado di verbosita', da 0 [default: 0]
    int verbose;
    ///Copia della riga di comando
    std::string command_line;
    //bool translate;
    ///Mostrare la demo di operazioni tra partizioni
    bool demo;

} options;

//start and load functions
void set_program_options(options &opt, int argc, char**argv);
void fill_entries_randomly(std::string *entries);
void fill_seq_from_file(options &opts, std::string* entries);
void generate_next_sequence(int *num_entry);
int load_config(options &opts, int *num_entry);
//more general 
///Coppia <entropia_shannon, numero_atomi>
typedef std::pair<double,int> entropy_pair;
template <typename T>  entropy_pair ordered_vector_entropy(const T *temp, int N);
template <typename partition_t> void calcola_matrice_distanze(const partition_t *X);
template <typename data_t> void print_array(const data_t *array, int len, const char *nome);
template <typename T> void ppmout(const T *grid, int sz, const char *filename);
template <typename T, typename U> void ppmout2(const T *grid1, const U* grid2, int sz, const char *filename);

#endif
