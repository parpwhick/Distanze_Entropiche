/*
 * File:   dopon_problem.cpp
 * Author: fake
 *
 * Created on September 30, 2012, 2:48 PM
 */

#include "dopon_problem.h"
#include <cmath>

#ifdef USE_EIGEN
template <typename spin_t> void dopon_problem<spin_t>::construct_problem_matrix(int spin) {
    H.resize(NN.N, NN.N);
    H.setZero();
    int z, sum_nn;

    for (int k = 0; k < NN.N; k++) {

        z = NN.fetch(k);
        sum_nn = 0;
        for (int m = 0; m < z; m++) {
            sum_nn += s[NN.vicini[m]];
            H(NN.vicini[m], k) = t;
        }

        H(k, k) = lambda * (s[k] == spin) + J * spin * 0.25 * sum_nn;

        if (confining==0)
                H(k, k) += -V * (2*(floor(k / L)/(L - 1))-1); //linear
        else
                H(k, k) += 3 * square(2*((floor(k / L)-confining) / (L - 1))); //parabola, centered in the middle row, y in [0,1]
    }
}
template void dopon_problem<char>::construct_problem_matrix(int spin);
template void dopon_problem<double>::construct_problem_matrix(int spin);

template <typename spin_t> double dopon_problem<spin_t>::calculate_lowest_energy(bool verbose) {
    SelfAdjointEigenSolver<MatrixXf> eigenproblem;

    double E_up = 1e4, E_down = 1e4;

    if (last_gs_spin != 1) {
        construct_problem_matrix(-1);
        eigenproblem.compute(H, EigenvaluesOnly);

        if (eigenproblem.info() == Success)
            E_down = eigenproblem.eigenvalues().coeff(0) / J;
        else
            fprintf(stderr, "Did not manage to find eigenvalue for spin down dopons\n");
    }

    if (last_gs_spin != -1) {
        construct_problem_matrix(+1);
        eigenproblem.compute(H, EigenvaluesOnly);
        if (eigenproblem.info() == Success)
            E_up = eigenproblem.eigenvalues().coeff(0) / J;
        else
            fprintf(stderr, "Did not manage to find eigenvalue for spin up dopons\n");
    }

    if (verbose)
        fprintf(stderr, "E(+) = %f, E(-) = %f\n", E_up, E_down);

    if (std::abs(E_down - E_up) > 1)
        last_gs_spin = (E_down < E_up) ? -1 : 1;
    else
        last_gs_spin = 0;

    last_energy = std::min(E_up, E_down);
    return last_energy;
}
template double dopon_problem<char>::calculate_lowest_energy(bool);
template double dopon_problem<double>::calculate_lowest_energy(bool);
#endif

//int vicini[4] = {(k - (k % L)+ ((k+L-1)%L)),((k/L)*L + ((k+L+1)%L)),(k+N-L)%N,(k+N+L)%N};
template <class spin_t> void dopon_problem<spin_t>::MultMv(double*in, double*out) {
    const int & spin = probed_spin;
    int z;
    double sum_nn;

    for (int k = 0; k < NN.N; k++) {
        //set the resulting vector element to 0
        out[k] = 0;
        //iterate over the neighbors
        z = NN.fetch(k);
        sum_nn = 0;
        for (int m = 0; m < z; m++) {
            sum_nn += s[NN.vicini[m]];
            //action of the "offdiagonal" terms, picking a t for every neighbor
            out[k] += in[NN.vicini[m]]; //*t;
        }
        //diagonal element, dopon spin dependent
        out[k] += (lambda * (s[k] == spin) + J * spin * 0.25 * sum_nn) * in[k];
        //voltage gradient
        if (confining==0)
                out[k] += -V * (floor(k / L)/(L - 1)) * in[k]; //linear
        else
                out[k] += 5 * square(2*((floor(k / L)-confining) / (L - 1))) * in[k]; //parabola, centered at "confining", y in [0,1]
    }
}
template void dopon_problem<char>::MultMv(double *in, double*out);
template void dopon_problem<double>::MultMv(double *in, double*out);

template <typename spin_t> double dopon_problem<spin_t>::lanczos_lowest_energy(bool verbose){
    double E_up=1e4, E_down=1e4;

    if (last_gs_spin != 1) {
        probed_spin = -1;
        E_down = lanczos_groundstate(false) / J;
    }

   if (last_gs_spin != -1) {
        probed_spin = +1;
        E_up = lanczos_groundstate(false) / J;
        
    }

    if (verbose)
        fprintf(stderr, "E(+) = %f, E(-) = %f\n", E_up, E_down);
    
    if(std::abs(E_down - E_up) > 1)
        last_gs_spin = (E_down < E_up) ? -1 : 1;
    else
        last_gs_spin = 0;

    last_energy = std::min(E_up,E_down);
    return last_energy;
}
template double dopon_problem<double>::lanczos_lowest_energy(bool);
template double dopon_problem<char>::lanczos_lowest_energy(bool);

///Normalize the input vector, v/sqrt(<v,v>) and return the norm <v,v>
double normalize(double *v, size_t L){
    double norm=0.0, norm2=0.0;
    for(size_t i=0; i<L;i++)
        norm2+=v[i]*v[i];
    norm = std::sqrt(norm2);
    for(size_t i=0; i<L;i++)
        v[i]/=norm;
    return norm;
}

double scalar_product(const double *v1, const double *v2, size_t L){
    double prod=0.0;
    for(size_t i=0; i<L;i++)
        prod+=v1[i]*v2[i];
    return prod;
}

void subtract_scaled3(double *v2,double a,const double *v1,double n, const double *v0,size_t L){
    for(size_t i=0; i<L;i++)
        v2[i]-= a * v1[i] + n * v0[i];
}

void subtract_scaled2(double *v2,double a,const double *v1,size_t L){
    for(size_t i=0; i<L;i++)
        v2[i]-= a * v1[i];
}

void sum_scaled(double *v2,double a,const double *v1,size_t N){
    for(size_t i=0; i<N;i++)
        v2[i]+= a * v1[i];
}

template <typename spin_t> double dopon_problem<spin_t>::lanczos_groundstate(int verbose) {
    int nconv;
    if(ground_state.empty())
        //we assign a uniform initial state, although it's a terrible choice
        //in certain cases the true groundstate could be the uniform vector,
        //leading to errors and nonconvergence in the Lanczos algorithm.
        ground_state.assign(NN.N, 1.0);
    if(a.empty()){
        a.resize(maxiterations);
        n.resize(maxiterations);
        d.resize(maxiterations);
        e.resize(maxiterations);
        eigenvectors.resize(maxiterations*maxiterations);
        fullvectors.resize(maxiterations * NN.N);
    }
    double *phi0;
    double *phi1;
    double *phi2;

    int N = NN.N; //vector length
    int m = 1; //iteration counter
    double gs, gs_old = 0.0;

    phi0 = &fullvectors[0];
    std::copy(ground_state.begin(), ground_state.end(),&phi0[0]);
    double gs_norm = normalize(phi0,N);
    if (verbose>1) printf("Initial state normalized with: %f\n", gs_norm);

    phi1 = &fullvectors[1*N];
    MultMv(phi0, phi1);
    a[0] = scalar_product(phi0, phi1,N);
    subtract_scaled2(phi1, a[0], phi0,N);
    n[0] = normalize(phi1,N);
    if (verbose>1) fprintf(stdout, "Iteration 0: a = %f, n = %f\n", a[0], n[0]);

    for (m = 1; m < maxiterations-1; m++) {
        phi2 = &fullvectors[(m+1)*N];
        phi1 = &fullvectors[m*N];
        phi0 = &fullvectors[(m-1)*N];
        
        MultMv(phi1, phi2);
        a[m] = scalar_product(phi2, phi1,N);
        subtract_scaled3(phi2, a[m], phi1, n[m - 1], phi0,N);
        n[m] = normalize(phi2,N);
        if (verbose>1) fprintf(stdout, "Iteration %d: a = %f, n = %f\n", m, a[m], n[m]);
        
        if (m > 18) {
            d.assign(&a[0], &a[m]);
            e.assign(&n[0], &n[m]);
            nconv = tqlrat(m, d.data(), e.data());
            if (nconv != 0)
                fprintf(stderr, "TQLRAT ERROR: bad convergence, diagonalized only %d values\n", nconv);
            gs = d[0];
            if (std::abs(gs - gs_old) < tol)//check for convergence           
                break;
            else
                gs_old = gs;
        }
    }
    
    if (m >= maxiterations-1) {
        fprintf(stderr, "LANCZOS: convergence not reached after %d steps, retrying from last known state\n", m);
        ground_state.assign(&phi2[0],&phi2[N]);
        return lanczos_groundstate(verbose);
    } else
        m--; //m is 1 higher than the last determined coefficient, needs to be -1
    
    if (verbose) {
        fprintf(stdout, "Mylanczos: convergence in %d iterations.\nE=[", m+1);
        for (int i = 0; i < std::min(m, 20); i++)
            fprintf(stdout, "%f, ", d[i]);
        fprintf(stdout, "]\n\n");
    }

    //now, the ground state eigenvector...
    
    //set the eigenvector matrix diagonal to 1
    //tql2 will reconstruct the eigenvector matrix
    eigenvectors.assign(eigenvectors.size(),0.0);
    for (int i = 0; i < m; i++)
        eigenvectors[i * m + i] = 1.0;
    d.assign(&a[0], &a[m]);
    e.assign(&n[0], &n[m]);
    //full diagonalization with eigenvectors
    nconv = tql2(m, d.data(), e.data(), eigenvectors.data());
    if (nconv != 0)
        fprintf(stderr, "TQL2 ERROR: bad convergence, diagonalized only %d vectors\n", nconv);
   
    ground_state.assign(N, 0.0);
    for (int j = 0; j < m; j++)
        //every component of the eigenvector in the Lanczos basis (eigenvectors[j])
        //multiplies the vectors forming said basis (phi1 at every iteration)
        sum_scaled(&ground_state[0], eigenvectors[j], &fullvectors[j*N], N);
    
    
    gs_norm = normalize(ground_state.data(),N);
    if (std::abs(gs_norm - 1) > 1e-7) {
        fprintf(stderr, "LANCZOS algorithm failure: for the groundstate after %d steps, |gs| = %f\n", m,gs_norm);
        fprintf(stderr, "Restart with a random initial state and check the Hamiltonian isn't constant\n");
        exit(1);
        
    }
    return d[0];
}
template double dopon_problem<double>::lanczos_groundstate(int verbose);
template double dopon_problem<char>::lanczos_groundstate(int verbose);


template <typename T> int signof(T number){
    return 2*(number >= 0)-1;
}

double pythag(double a, double b){
//      PYTHAG = sqrt ( A * A + B * B )
    double p,r,s,t,u;
    p = std::max(std::abs(a), std::abs(b));

    if (p != 0.0) {
        r = std::min(std::abs(a), std::abs(b)) / p;
        r = r * r;

        while (1) {
            t = 4.0 + r;
            if (t == 4.0)
                break;
            s = r / t;
            u = 1.0 + 2.0 * s;
            p = u * p;
            r = (s / u) * (s / u) * r;
        }
    }
    return p;
}

int tql2(int n, double d[], double e[], double z[]){
//****************************************************************************80
//
//  Purpose:
//
//    TQL2 computes all eigenvalues/vectors, real symmetric tridiagonal matrix.
//
//  Discussion:
//
//    This subroutine finds the eigenvalues and eigenvectors of a symmetric
//    tridiagonal matrix by the QL method.  The eigenvectors of a full
//    symmetric matrix can also be found if TRED2 has been used to reduce this
//    full matrix to tridiagonal form.
//
//  Author:
//
//    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
//    Klema, Moler.
//    C++ version by John Burkardt.
//

//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input/output, double D[N].  On input, the diagonal elements of
//    the matrix.  On output, the eigenvalues in ascending order.  If an error
//    exit is made, the eigenvalues are correct but unordered for indices
//    1,2,...,IERR-1.
//
//    Input/output, double E[N].  On input, E(1:N-1) contains the
//    subdiagonal elements of the input matrix, and E(N) is arbitrary.
//    On output, E has been destroyed.
//
//    Input, double Z[N*N].  On input, the transformation matrix
//    produced in the reduction by TRED2, if performed.  If the eigenvectors of
//    the tridiagonal matrix are desired, Z must contain the identity matrix.
//    On output, Z contains the orthonormal eigenvectors of the symmetric
//    tridiagonal (or full) matrix.  If an error exit is made, Z contains
//    the eigenvectors associated with the stored eigenvalues.
//
//    Output, int TQL2, error flag.
//    0, normal return,
//    J, if the J-th eigenvalue has not been determined after
//    30 iterations.
//
    double c=0,c2=0,c3=0;
    double dl1, el1, f, g, h;
    int i, ierr, ii, j, k, l, l1, l2, m, mml;
    double p, r, s, s2=0, t;
    double tst1, tst2;

    ierr = 0;
    if (n == 1) 
        return ierr;

    f = 0.0;
    tst1 = 0.0;
    e[n - 1] = 0.0;

    for (l = 0; l < n; l++) {
        j = 0;
        h = std::abs(d[l]) + std::abs(e[l]);
        tst1 = std::max(tst1, h);
        //
        //  Look for a small sub-diagonal element.
        //
        for (m = l; m < n; m++) {
            tst2 = tst1 + std::abs(e[m]);
            if (tst2 == tst1) {
                break;
            }
        }

        if (m != l) {
            for (;;) {
                if (30 <= j) {
                    ierr = l + 1;
                    return ierr;
                }

                j = j + 1;
                //
                //  Form shift.
                //
                l1 = l + 1;
                l2 = l1 + 1;
                g = d[l];
                p = (d[l1] - g) / (2.0 * e[l]);
                r = pythag(p, 1.0);
                d[l] = e[l] / (p + signof(p) * std::abs(r));
                d[l1] = e[l] * (p + signof(p) * std::abs(r));
                dl1 = d[l1];
                h = g - d[l];
                for (i = l2; i < n; i++) {
                    d[i] = d[i] - h;
                }
                f = f + h;
                //
                //  QL transformation.
                //
                p = d[m];
                c = 1.0;
                c2 = c;
                el1 = e[l1];
                s = 0.0;
                mml = m - l;

                for (ii = 1; ii <= mml; ii++) {
                    c3 = c2;
                    c2 = c;
                    s2 = s;
                    i = m - ii;
                    g = c * e[i];
                    h = c * p;
                    r = pythag(p, e[i]);
                    e[i + 1] = s * r;
                    s = e[i] / r;
                    c = p / r;
                    p = c * d[i] - s * g;
                    d[i + 1] = h + s * (c * g + s * d[i]);
                    //
                    //  Form vector.
                    //
                    for (k = 0; k < n; k++) {
                        h = z[k + (i + 1) * n];
                        z[k + (i + 1) * n] = s * z[k + i * n] + c * h;
                        z[k + i * n] = c * z[k + i * n] - s * h;
                    }
                }
                p = -s * s2 * c3 * el1 * e[l] / dl1;
                e[l] = s * p;
                d[l] = c * p;
                tst2 = tst1 + std::abs(e[l]);

                if (tst2 <= tst1) 
                    break;
            }
        }
        d[l] = d[l] + f;
    }
    //
    //  Order eigenvalues and eigenvectors.
    //
    for (ii = 1; ii < n; ii++) {
        i = ii - 1;
        k = i;
        p = d[i];
        for (j = ii; j < n; j++) {
            if (d[j] < p) {
                k = j;
                p = d[j];
            }
        }

        if (k != i) {
            d[k] = d[i];
            d[i] = p;
            for (j = 0; j < n; j++) {
                t = z[j + i * n];
                z[j + i * n] = z[j + k * n];
                z[j + k * n] = t;
            }
        }
    }
    return ierr;
}
//****************************************************************************80

int tqlrat(int n, double d[], double e[]){
//****************************************************************************80
//  Purpose:
//
//    TQLRAT computes all eigenvalues of a real symmetric tridiagonal matrix.
//
//  Discussion:
//
//    This subroutine finds the eigenvalues of a symmetric
//    tridiagonal matrix by the rational QL method.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input/output, double D[N].  On input, D contains the diagonal
//    elements of the matrix.  On output, D contains the eigenvalues in ascending
//    order.  If an error exit was made, then the eigenvalues are correct
//    in positions 1 through IERR-1, but may not be the smallest eigenvalues.
//
//    Input/output, double E[N], contains in positions 0 through N-2
//    the subdiagonal elements of the matrix.  E(N-1) is
//    arbitrary.  On output, E has been overwritten by workspace
//    information.
//
//    Output, int TQLRAT, error flag.
//    0, for no error,
//    J, if the J-th eigenvalue could not be determined after 30 iterations.
    double b=0, c=0,f,g,h;
    int i, ierr,ii,j,l,l1,m=0,mml;
    double p,r,s,t;

    ierr = 0;
    if (n == 1) 
        return ierr;

    f = 0.0;
    t = 0.0;
    for(i=0;i<n;i++)
        e[i] *= e[i];
    e[n - 1] = 0.0;

    for (l = 0; l < n; l++) {
        j = 0;
        h = std::abs(d[l]) + sqrt(e[l]);

        if (t <= h) {
            t = h;
            b = std::abs(t) * epsilon;
            c = b * b;
        }
        //
        //  Look for small squared sub-diagonal element.
        //
        for (m = l; m < n; m++)
            if (e[m] <= c)
                break;
        
        if (m != l) {
            for (;;) {
                if (30 <= j) {
                    ierr = l + 1;
                    return ierr;
                }

                j = j + 1;
                //
                //  Form shift.
                //
                l1 = l + 1;
                s = sqrt(e[l]);
                g = d[l];
                p = (d[l1] - g) / (2.0 * s);
                r = pythag(p, 1.0);
                d[l] = s / (p + std::abs(r) * signof(p));
                h = g - d[l];
                for (i = l1; i < n; i++) {
                    d[i] = d[i] - h;
                }
                f = f + h;
                //
                //  Rational QL transformation.
                //
                g = d[m];
                if (g == 0.0)
                    g = b;

                h = g;
                s = 0.0;
                mml = m - l;

                for (ii = 1; ii <= mml; ii++) {
                    i = m - ii;
                    p = g * h;
                    r = p + e[i];
                    e[i + 1] = s * r;
                    s = e[i] / r;
                    d[i + 1] = h + s * (h + d[i]);
                    g = d[i] - e[i] / g;
                    if (g == 0.0)
                        g = b;
                    h = g * p / r;
                }
                e[l] = s * g;
                d[l] = h;
                //
                //  Guard against underflow in convergence test.
                //
                if (h == 0.0)
                    break;

                if (std::abs(e[l]) <= std::abs(c / h))
                    break;

                e[l] = h * e[l];

                if (e[l] == 0.0)
                    break;
            }
        }

        p = d[l] + f;
        //
        //  Order the eigenvalues.
        //
        for (i = l; 0 <= i; i--) {
            if (i == 0) {
                d[i] = p;
                break;
            } else if (d[i - 1] <= p) {
                d[i] = p;
                break;
            }
            d[i] = d[i - 1];
        }
    }

    return ierr;
}