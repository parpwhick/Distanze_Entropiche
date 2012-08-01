/* 
 * File:   strutture.h
 * Author: fake
 *
 * Created on January 13, 2012, 10:51 PM
 */

#ifndef STRUTTURE_H
#define	STRUTTURE_H

#include <vector>
#include <string>
#include <stdint.h>

typedef int32_t label_t;
typedef uint64_t product_t;

enum simulation_t {
    METROPOLIS,
    MICROCANONICAL
};

enum part_type {
    GENERAL_PARTITION,
    LINEAR_PARTITION
};

enum DISTANZE {
    SHAN = 1,
    TOP = 1 << 1,
    RID = 1 << 2,
    RID_TOP = 1 << 3,
    HAMM = 1 << 20
};

enum source {
    LINEARE,
    RETICOLO_2D,
    SIERPINSKI,
    FUZZY,
    RANDOM,
    FROM_FILE,
    SIMULATION,
};

enum red_strategy{
    COMUNE,
    DIRETTA
};

typedef struct {
    int seq_len;
    int n_seq;
    int lato;
    int sierpinski_gen;
    int n_symbols;
    part_type partition_type;
    int epsilon;

    source letto_da;
    source topologia;
    char state_filename[255];
    char adj_vec_1[255];
    char adj_vec_2[255];
    int fuzzy;
    red_strategy riduzione;

    int da_calcolare;
    bool write;
    bool distance;
    int threads;
    bool graphics;

    int seed;
    simulation_t simulation_type;
    int steps;
    double beta;
    int max_link_energy;
    
    int verbose;
    //bool translate;

} options;

class general_partition;
class linear_partition;


//start and load functions
void set_program_options(options &opt, int argc, char**argv);
void fill_entries_randomly(std::string *entries);
void fill_seq_from_file(options &opts, std::string* entries);
void generate_next_sequence(int *num_entry);
int load_config(options &opts, int *num_entry);

//more general 
typedef std::pair<double,int> entropy_pair;
template <typename T>  entropy_pair ordered_vector_entropy(const T *temp, int N);
template <typename partition_t> void calcola_matrice_distanze(const partition_t *X);
template <typename data_t> void print_array(const data_t *array, int len, const char *nome);
template <typename T> void ppmout(const T *grid, int sz, const char *filename);
template <typename T, typename U> void ppmout2(const T *grid1, const U* grid2, int sz, const char *filename);

#endif
