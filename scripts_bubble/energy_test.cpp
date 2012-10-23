#include <iostream>
#include <fstream>
#include <iterator>
#include <vector>
#include <algorithm>
#include <cmath>

#include "dopon_problem.h"

using namespace std;
int main(int argc, char**argv){
 /*   char *where;
    if (argc>1)
        where = argv[1];
    else
        where = "spins.txt";

  std::ifstream in(where);

  if (!in) {
    std::cout << "Cannot open file.\n";
    return 0;
  }
*/
    typedef double input_t;
  std::istream_iterator<input_t> ii(cin);
  std::istream_iterator<input_t> eos;

  std::vector<input_t> spins(ii,eos);

  int N = spins.size();
  int L = static_cast<int>(std::round(std::sqrt(N)));
  if(L*L != N){
      cerr << "wrong matrix size (i.e. not square)" << endl;
      return 1;
  }
  cout << "given a " << L << "x" << L << " matrix" << endl << endl;

  adj_struct top = adiacenza_toroidal_lattice(L);
  dopon_problem<input_t> prob(1,1,300,top);
  prob.set_spin_array(spins.data());

  cout << "Lanczos:" << endl;
  prob.lanczos_lowest_energy(true);
  cout << endl;
  cout << "Full:" << endl;
  prob.calculate_lowest_energy(true);
}
