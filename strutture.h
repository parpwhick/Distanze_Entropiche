/* 
 * File:   strutture.h
 * Author: fake
 *
 * Created on January 13, 2012, 10:51 PM
 */

#ifndef STRUTTURE_H
#define	STRUTTURE_H

#include "rand55.h"
#include "adj_handler.h"
#include <string>
#include <stdint.h>

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
    FUZZY,
    RANDOM,
    FROM_FILE,
    SIMULATION,
};

typedef struct {
    int seq_len;
    int n_seq;
    int lato;
    int n_symbols;
    part_type partition_type;

    source letto_da;
    source topologia;
    char state_filename[255];
    char adj_vec_1[255];
    char adj_vec_2[255];
    int fuzzy;

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

typedef struct {
    int value;
    int pos;
} labelled_array;

typedef int32_t label_t;

class basic_partition {
public:
    //number of atoms found
    label_t n;
    //total length of the partition
    label_t N;

    double entropia_shannon;
    double entropia_topologica;
};

class linear_partition : public basic_partition {
public:
    //binary atom beginning codes {1,0,0,0,1,1,...}
    int *binary;

    template <typename T> linear_partition(const T*seq, int len);
    template <typename T> void fill(const T*seq, int len);
    void print();

    linear_partition(int len = 0) {
        n = 0;
        N = len;
        binary = new int[len];
    }

    ~linear_partition() {
        if (N)
            delete[] binary;
        n = N = 0;
    }
};

class atom {
public:
    label_t size;
    label_t start;
    label_t end;

    atom() {
        size = 0;
        end = 0;
        start = 0;
    }
};

class general_partition : public basic_partition {
private:
    //nr. di atomi allocati per questa partizione (ottimizzazione)
    label_t allocated_n;

    void allocate(label_t len);
    void allocate_atoms(label_t n1);
    void sort_entropy();
    void relabel();
    //nearest neighbors
    label_t *prev_site;

    atom * atomi;
public:
    //labels identify generic atoms across the partition
    label_t *labels;

    void from_nnb(label_t **NNB);
    void from_configuration(const int *configuration, adj_struct adj, int N1=0);

    general_partition(int len = 0);
    ~general_partition();
    void reduce(const general_partition &ridurre, const general_partition &common);
    void linear_intersection(const general_partition &p1, const general_partition &p2);

    void print();
    void print_cluster_adjacency();

    atom& find_atom(const atom &atomo1) const {
        return atomi[labels[atomo1.start]];
    }

    class Iterator {
    private:
        label_t _site;
        const label_t *_next;

    public:

        Iterator(int dove, const label_t *vicini) : _site(dove), _next(vicini) {
        };

        int operator*() {
            return (_site);
        }

        bool operator==(Iterator due) {
            return (_site == due._site);
        }

        bool operator!=(Iterator due) {
            return (_site != due._site);
        }

        int operator++() {
            if (_site == _next[_site])
                _site = -1;
            else
                _site = _next[_site];
            return (_site);
        }

        int operator++(int) {
            return (operator++());
        }

    };

    Iterator begin(const int where) const {
        return Iterator(atomi[where].end, prev_site);
    }

    Iterator begin(const atom &where) const {
        return Iterator(where.end, prev_site);
    }

    Iterator end() const {
        return Iterator(-1, prev_site);
    }
};

typedef general_partition::Iterator Iter_t;

class distance {
private:
    void allocate(int n);

    int *common_factor;
    int *reduced1;
    int *reduced2;
    int *product_reduced;
    uint64_t *product;
    general_partition ridotto1;
    general_partition ridotto2;
    general_partition partizione_comune;

    int N;

public:
    double dist_shan;
    double dist_shan_r;
    double dist_top;
    double dist_top_r;
    double dist_ham;

    void dist(const linear_partition &p1, const linear_partition &p2);
    void dist(const general_partition &p1, const general_partition &p2);
    void fill(const linear_partition& e1, const linear_partition& e2);
    void fill(const general_partition& e1, const general_partition& e2);
    template <typename T> void hamming_distance(const T* seq1, const T* seq2);

    ~distance();
    distance(int n = 1000);
    distance(const distance &d1);


};
//start and load functions
void set_program_options(options &opt, int argc, char**argv);
void fill_entries_randomly(std::string *entries);
void fill_seq_from_file(options &opts, std::string* entries);
void generate_next_sequence(int *num_entry);
int load_config(options &opts, int *num_entry);

//more general 
template <typename T> double entropy_binary_partition(const T *p, int n);
template <typename partition_t> void calcola_matrice_distanze(const partition_t *X);
template <typename data_t> void print_array(const data_t *array, int len, const char *nome);
template <typename T> void ppmout(const T *grid, int sz, const char *filename);
template <typename T, typename U> void ppmout2(const T *grid1, const U* grid2, int sz, const char *filename);

#endif