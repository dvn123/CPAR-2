#include <iostream>
#include <limits>
#include <sstream>
#include <cmath>
#include <vector>
#include <omp.h>

using namespace std;

bool debug_msg_i = false;
bool debug_msg_j = false;
bool print_primes = false;

void sieve(long limit, vector<bool> &is_prime) {
  long sq = sqrt(limit*2);
  //#pragma omp parallel for shared(is_prime) schedule(dynamic)
  for(long i = 1; i <= (sq-1)/2; i++) {
    debug_msg_i && cout << "i - " << i << endl;
    debug_msg_i && cout << 1+2*i << " < " << sqrt(limit*2) << endl;
    if(is_prime[i]) {
      debug_msg_i && cout << "Chosen" << endl;
      for (long j = i + (i*(2*i+1)); j <= limit; j += 2*i+1) {
        debug_msg_j && cout << "j - " << j << endl;
        is_prime[j] = false;
      }
    }
  }
}

int main(int argc, char *argv[]) {
  long limit = -1;
  if (argc == 2) {
    stringstream ss(argv[--argc]);
    ss >> limit;

    if (limit < 1 || ss.fail()) {
      cerr << "USAGE:\n  sieve LIMIT\nwhere LIMIT in the range [1, " << numeric_limits<long>::max() << "]" << endl;
      return 2;
    }
  } else {
    cerr << "USAGE:\n  sieve LIMIT\n\nwhere LIMIT in the range [1, " << numeric_limits<long>::max() << "]" << endl;
    return 2;
  }

  //omp_set_num_threads(omp_get_num_procs());

  long limit_t = floor(limit/2.0);

#pragma omp parallel
  printf("Hello, world.\n");


  vector<bool> is_prime(limit_t, true);

  clock_t begin_time = clock();
  sieve(limit_t, is_prime);

  cout << "Time: " << float(clock () - begin_time)/CLOCKS_PER_SEC;

  //Print
  print_primes && cout << endl << "Primes up to " << limit << ":" << endl;
  for(int i = 0; i < is_prime.size(); i++) {
    print_primes && is_prime[i] && cout << 2*i+1 << " ";
  }
}