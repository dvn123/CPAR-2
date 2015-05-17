#include <iostream>
#include <limits>
#include <sstream>
#include <cmath>
#include <vector>

#define REAL_NUMBER(n) \
(1+2*n)

#define ARRAY_INDEX(n) \
((n-1)/2)

using namespace std;

bool debug_msg_i = false;
bool debug_msg_j = false;

void sieve(int limit, vector<bool> &is_prime) {
  for(int i = 1; REAL_NUMBER(i)*REAL_NUMBER(i) <= limit*2; i++) {
    debug_msg_i && cout << "i - " << i << endl;
    debug_msg_i && cout << REAL_NUMBER(i) << " < " << sqrt(limit*2) << endl;
    if(is_prime[i]) {
      debug_msg_i && cout << "Chosen" << endl;
      for (int j = i + (i*(REAL_NUMBER(i))); j <= limit; j += REAL_NUMBER(i)) {
        debug_msg_j && cout << "j - " << j << endl;
        is_prime[j] = false;
      }
    }
  }
}

int main(int argc, char *argv[]) {
  int limit = -1;
  if (argc == 2) {
    stringstream ss(argv[--argc]);
    ss >> limit;

    if (limit < 1 || ss.fail()) {
      cerr << "USAGE:\n  sieve LIMIT\nwhere LIMIT in the range [1, " << numeric_limits<int>::max() << "]" << endl;
      return EXIT_FAILURE;
    }
  } else {
    cerr << "USAGE:\n  sieve LIMIT\n\nwhere LIMIT in the range [1, " << numeric_limits<int>::max() << "]" << endl;
    return EXIT_FAILURE;
  }

  int limit_t = floor(limit/2.0);

  vector<bool> is_prime(limit_t, true);

  sieve(limit_t, is_prime);

  //Print
  cout << "Primes up to " << limit << ":" << endl;
  for(int i = 0; i < is_prime.size(); i++) {
    is_prime[i] && cout << 2*i+1 << " ";
    //cout << is_prime[i] << endl;
  }

  return EXIT_SUCCESS;
}