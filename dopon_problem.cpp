/*
 * File:   dopon_problem.cpp
 * Author: fake
 *
 * Created on September 30, 2012, 2:48 PM
 */

#include "dopon_problem.h"

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

    if (abs(E_down - E_up) > 1)
        last_gs_spin = (E_down < E_up) ? -1 : 1;
    else
        last_gs_spin = 0;

    last_energy = std::min(E_up, E_down);
    return last_energy;
}
template double dopon_problem<char>::calculate_lowest_energy(bool);
template double dopon_problem<double>::calculate_lowest_energy(bool);
#endif USE_EIGEN

template <class spin_t> template <typename T> void dopon_problem<spin_t>::MultMv(T*in, T*out) {

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
            out[k] += in[NN.vicini[m]] * t;
        }
        //diagonal element, dopon spin dependent
        out[k] += (lambda * (s[k] == spin) + J * spin * 0.25 * sum_nn) * in[k];
        //voltage gradient
        if (confining==0)
                out[k] += -V * (2*(floor(k / L)/(L - 1))-1) * in[k]; //linear
        else
                out[k] += 5 * square(2*((floor(k / L)-confining) / (L - 1))) * in[k]; //parabola, centered at "confining", y in [0,1]
    }
}
template void dopon_problem<char>::MultMv(double *in, double*out);
template void dopon_problem<double>::MultMv(double *in, double*out);

template<class FLOAT, class EIGPROB>
void Solution(dopon_problem<char> &A, EIGPROB &Prob)
/*
  This function prints eigenvalues and eigenvetors on standard "cout"
  stream and exemplifies how to retrieve information from ARPACK++ classes.
*/

{
    using namespace std;
  int   i, n, nconv, mode;
  FLOAT *Ax;
  FLOAT *ResNorm;

  /*
     ARPACK++ includes some functions that provide information
     about the problem. For example, GetN furnishes the dimension
     of the problem and ConvergedEigenvalues the number of
     eigenvalues that attained the required accuracy. GetMode
     indicates if the problem was solved in regular,
     shift-and-invert or other mode.
  */

  n     = Prob.GetN();
  nconv = Prob.ConvergedEigenvalues();
  mode  = Prob.GetMode();
/*
  cout << endl << endl << "Testing ARPACK++ class ARSymEig \n";
  cout << "Real symmetric eigenvalue problem: A*x - lambda*x" << endl;
  switch (mode) {
  case 1:
    cout << "Regular mode" << endl << endl;
    break;
  case 3:
    cout << "Shift and invert mode" << endl << endl;
  }

  cout << "Dimension of the system            : " << n             << endl;
  cout << "Number of 'requested' eigenvalues  : " << Prob.GetNev() << endl;
  cout << "Number of 'converged' eigenvalues  : " << nconv         << endl;
  cout << "Number of Arnoldi vectors generated: " << Prob.GetNcv() << endl;
  cout << endl;
*/
  /*
    EigenvaluesFound is a boolean function that indicates
    if the eigenvalues were found or not. Eigenvalue can be
    used to obtain one of the "converged" eigenvalues. There
    are other functions that return eigenvectors elements,
    Schur vectors elements, residual vector elements, etc.
  */

  if (Prob.EigenvaluesFound()) {
    cout << "Eigenvalues:" << endl;
    for (i=0; i<nconv; i++) {
      cout << "  lambda[" << (i+1) << "]: " << Prob.Eigenvalue(i) << endl;
    }
    cout << endl;
  }

  /*
    EigenvectorsFound indicates if the eigenvectors are
    available. RawEigenvector is one of the functions that
    provide raw access to ARPACK++ output data. Other functions
    of this type include RawEigenvalues, RawEigenvectors,
    RawSchurVector, RawResidualVector, etc.
  */

  if (Prob.EigenvectorsFound()) {

    // Printing the residual norm || A*x - lambda*x ||
    // for the nconv accurately computed eigenvectors.

    Ax      = new FLOAT[n];
    ResNorm = new FLOAT[nconv];

    for (i=0; i<nconv; i++) {
      A.MultMv(Prob.RawEigenvector(i),Ax);
      axpy(n, -Prob.Eigenvalue(i), Prob.RawEigenvector(i), 1, Ax, 1);
      ResNorm[i] = nrm2(n, Ax, 1) / fabs(Prob.Eigenvalue(i));
    }

    for (i=0; i<nconv; i++) {
      cout << "||A*x(" << (i+1) << ") - lambda(" << (i+1);
      cout << ")*x(" << (i+1) << ")||: " << ResNorm[i] << "\n";
    }
    cout << "\n";

    delete[] Ax;
    delete[] ResNorm;

  }

} // Solution

template <typename spin_t> double dopon_problem<spin_t>::lanczos_lowest_energy(bool verbose){
    int nconv;
    double E_up=1e4, E_down=1e4;

    if(ground_state.empty())
        ground_state.assign(NN.N,1.0);

    static double tol = 1.0e-7;
    static int maxiterations = 5000;
    static int arnoldi_vectors = 7;
    ARSymStdEig<double, dopon_problem<spin_t> >
            dprob(NN.N, 1, this, &dopon_problem<spin_t>::MultMv, "SA", arnoldi_vectors, tol, maxiterations, ground_state.data());

    dprob.ChangeTol(1.0e-7);
    dprob.ChangeMaxit(5000);

    if (last_gs_spin != 1) {
        dprob.ChangeNcv(arnoldi_vectors);
        probed_spin = -1;
        nconv = dprob.FindEigenvectors();

        if (nconv)
            E_down = dprob.Eigenvalue(0) / J;
        else{
            fprintf(stderr, "Did not manage to find eigenvalue for spin down dopons\n");
            arnoldi_vectors += 3;
            return lanczos_lowest_energy();
        }
        ground_state.assign(dprob.RawEigenvector(0),dprob.RawEigenvector(0)+NN.N);
    }

    if (last_gs_spin != -1) {
        probed_spin = +1;
        dprob.ChangeNcv(arnoldi_vectors);
        nconv = dprob.FindEigenvectors();

        if (nconv)
            E_up = dprob.Eigenvalue(0) / J;
        else{
            fprintf(stderr, "Did not manage to find eigenvalue for spin up dopons\n");
            arnoldi_vectors += 3;
            return lanczos_lowest_energy();
        }
        ground_state.assign(dprob.RawEigenvector(0),dprob.RawEigenvector(0)+NN.N);
    }

    if (verbose)
        fprintf(stderr, "E(+) = %f, E(-) = %f\n", E_up, E_down);

    if(abs(E_down - E_up) > 1)
        last_gs_spin = (E_down < E_up) ? -1 : 1;
    else
        last_gs_spin = 0;

    last_energy = std::min(E_up,E_down);
    return last_energy;
}
template double dopon_problem<char>::lanczos_lowest_energy(bool);
template double dopon_problem<double>::lanczos_lowest_energy(bool);