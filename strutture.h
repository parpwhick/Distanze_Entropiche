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
#include <vector>
#include <string>
#include <iostream>

#include <algorithm>
#include <stdio.h>

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
    SIERPINSKI,
    FUZZY,
    RANDOM,
    FROM_FILE,
    SIMULATION,
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
};

class general_partition : public basic_partition {
private:
    //nr. di atomi allocati per questa partizione (ottimizzazione)
    label_t allocated_n;

    void allocate(label_t len);
    void entropy_calculation();
    void relabel();
    //nearest neighbors
    std::vector<label_t> prev_site;
    std::vector<atom> atomi;
    general_partition(const general_partition &);
    general_partition & operator=(const general_partition &);
public:
    //labels identify generic atoms across the partition
    std::vector<label_t> labels;

    general_partition(int len = 0);

    template <typename data_t> void from_configuration(const data_t *configuration, const adj_struct & adj, int N1=0);
    void reduce(const general_partition &ridurre, const general_partition &common);
    void linear_intersection(const general_partition &p1, const general_partition &p2);
    void product(const general_partition & p1, const general_partition & p2);

    void print();
    void print_cluster_adjacency();

    const atom& find_atom(const atom &atomo1) const {
        return atomi[labels[atomo1.start]];
    }
    const atom& find_atom(const label_t inizio) const {
        return atomi[labels[inizio]];
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
        return Iterator(atomi[where].end, & prev_site[0]);
    }

    Iterator begin(const atom &where) const {
        return Iterator(where.end, & prev_site[0]);
    }

    Iterator end() const {
        return Iterator(-1, 0);
    }
};

typedef general_partition::Iterator Iter_t;

class distance {
private:
    void allocate(int n);

    std::vector<int> common_factor;
    std::vector<int> reduced1;
    std::vector<int> reduced2;
    std::vector<int> product_reduced;
    std::vector<uint64_t> product;
    general_partition ridotto1;
    general_partition ridotto2;
    general_partition partizione_comune;

    int N;
    void calc_distance(const general_partition &p1, const general_partition &p2);

public:
    double dist_shan;
    double dist_shan_r;
    double dist_top;
    double dist_top_r;
    double dist_ham;

    void dist(const linear_partition& e1, const linear_partition& e2);
    void dist(const general_partition& e1, const general_partition& e2);
    void operator()(const general_partition& e1, const general_partition& e2){
        dist(e1,e2);
    }
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