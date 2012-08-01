/* 
 * File:   partizioni.h
 * Author: fake
 *
 * Created on August 3, 2012, 2:39 PM
 */

#ifndef PARTIZIONI_H
#define	PARTIZIONI_H

#include <vector>
#include <string>
#include <stdint.h>
#include "adj_handler.h"

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

    template <typename data_t> void from_configuration(const data_t *configuration, const adj_struct & adj);
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

    class Iterator : public std::iterator_traits<label_t>{

    private:
        label_t _site;
        const label_t *_next;

    public:
        typedef label_t value_type;
        typedef label_t *pointer;
        typedef label_t &reference;
        typedef label_t difference_type;
    typedef std::forward_iterator_tag iterator_category;
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



#endif	/* PARTIZIONI_H */

