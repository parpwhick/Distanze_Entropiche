/* 
 * File:   dopon_problem.h
 * Author: fake
 *
 * Created on September 30, 2012, 2:48 PM
 */

#ifndef DOPON_PROBLEM_H
#define	DOPON_PROBLEM_H

#include "adj_handler.h"
#include "blas1c.h"
#include "arssym.h"
#include <iostream>
#include <cstdio>

#ifdef USE_EIGEN
#include "eigen3/Eigen/Dense"
using namespace Eigen;
#endif

template <typename T> T square (T val){
    return val*val;
}

template <typename spin_t>
class dopon_problem {
public:

    dopon_problem(double t1, double J1, double lambda1, const adj_struct & adj) :
    NN(adj), t(t1), J(J1),  lambda(lambda1)  {
        last_energy = 0;
        last_gs_spin = 0;
        confining = 0;
        V = 0;
    };

    double calculate_lowest_energy(bool verbose = false);
    double lanczos_lowest_energy(bool verbose=false);
    //double* get_ground_state();
    
    void set_spin_array(const spin_t *spin_array) {
        s = spin_array;
    }

    int lowest_energy_spin() {
        return last_gs_spin;
    }

    void set_t(double t_new) {
        t = t_new;
    }

    void set_lambda(double lambda_new) {
        lambda = lambda_new;
    }

    void set_J(double J_new){
        J = J_new;
    }
    
    void set_V(double V_new){
        V = V_new;
        lanczos_lowest_energy();
    }
    
    void set_confining(int P_new){
        confining = P_new;
        lanczos_lowest_energy();
    }
    
    void set_L(double L_new){
        L = L_new;
    }

    double last_energy;
    vector<double> ground_state;

    template <typename T> void MultMv (T*in, T*out);
    template <typename T> void nutrlan_multiply(const int vector_length, const int n_vectors, const double *in, const int in_vec_size,
	     double *out, const int out_vec_size, void*);
private:
    const adj_struct & NN;
#ifdef USE_EIGEN
    MatrixXf H;
#endif

    const spin_t *s;
    double t;
    double J;
    double lambda;
    
    int L;
    int last_gs_spin;
    int probed_spin;    

    double V;
    int confining;
    void construct_problem_matrix(int spin);
};

#endif	/* DOPON_PROBLEM_H */

