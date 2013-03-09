/* 
 * File:   dopon_problem.h
 * Author: fake
 *
 * Created on September 30, 2012, 2:48 PM
 */

#ifndef DOPON_PROBLEM_H
#define	DOPON_PROBLEM_H

#include "adj_handler.h"
#include <iostream>
#include <cstdio>
#include <vector>

#ifdef USE_EIGEN
#include "eigen3/Eigen/Dense"
using namespace Eigen;
#endif

template <typename T> T square (T val){
    return val*val;
}

int tqlrat(int n, double d[], double e2[]);
int tql2(int n, double d[], double e[], double z[]);
double pythag(double a, double b);

const double epsilon=2.220446049250313E-016;
const double lambda=300.0;
const double tol=1.0e-7;
const double maxiterations = 200;

template <typename spin_t>
class dopon_problem {
public:

    const adj_struct & NN;
    const spin_t *s;
    int L;
    int confining;
    double J;
    double V;
    double last_energy;
    int last_gs_spin;
    int probed_spin;
    
    
    vector<double> ground_state;
    vector<double> a; //Diagonal of tridiagonal matrix
    vector<double> d; //Eigenvalues in output
    vector<double> n; //Offdiagonal of tridiagonal matrix
    vector<double> e; //Working array for diagonalization
    vector<double> eigenvectors, fullvectors;
    
    dopon_problem(double J1, const adj_struct & adj) :
    NN(adj), J(J1) {
        last_energy = 0;
        last_gs_spin = 0;
        confining = 0;
        V = 0;
    };

    double calculate_lowest_energy(bool verbose = false);
    double lanczos_lowest_energy(bool verbose=false);
    void MultMv (double*in, double*out);
    double lanczos_groundstate(int verbose=0);
    
    void set_spin_array(const spin_t *spin_array) {
        s = spin_array;
    }

    int lowest_energy_spin() {
        return last_gs_spin;
    }

    void set_J(double J_new){
        J = J_new;
    }
    
    void set_V(double V_new){
        V = V_new;
        lanczos_lowest_energy();
    }
    
    void set_confining(double P_new){
        confining = P_new;
        lanczos_lowest_energy();
    }
    
    void set_L(double L_new){
        L = L_new;
    }
    
#ifdef USE_EIGEN
    MatrixXf H;
#endif
    void construct_problem_matrix(int spin);
};

#endif	/* DOPON_PROBLEM_H */

