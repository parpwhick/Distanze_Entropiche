/* 
 * File:   distance.h
 * Author: fake
 *
 * Created on August 3, 2012, 2:46 PM
 */

#ifndef DISTANCE_H
#define	DISTANCE_H

#include <vector>
#include <algorithm>
#include "partizioni.h"

class distance {
private:
    void allocate(int n);

    std::vector<char> common_factor;
    std::vector<char> reduced1;
    std::vector<char> reduced2;
    std::vector<char> product_reduced;
    std::vector<char> binary_product;
    std::vector<product_t> product;
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

#endif	/* DISTANCE_H */

