#include <iostream>
#include <limits>
#include <sstream>
#include <cmath>
#include <vector>
#include <omp.h>

#define BLOCK_LOW(id,p,n) \
((id)*(n)/(p))

#define BLOCK_HIGH(id,p,n) \
(BLOCK_LOW((id)+1,p,n)-1)

#define BLOCK_SIZE(id,p,n) \
(BLOCK_LOW((id)+1)-BLOCK_LOW(id))

#define BLOCK_OWNER(index,p,n) \
(((p)*(index)+1)-1)/(n))

using namespace std;

bool debug_msg_primes_used = true;
bool debug_msg_marking = true;
bool debug_msg_block_assignment = true;
bool debug_msg_first_index = true;

int n_threads;

void sieveBlockwise(int limit, vector<bool> &is_prime) {
  int prime_i = 1;
  #pragma omg parallel for default(shared) private(i)
  for(int i = 0; i <= n_threads - 1; i++) {
    int start = BLOCK_LOW(i, n_threads, limit);
    int end = BLOCK_HIGH(i,n_threads,limit);

    if(debug_msg_block_assignment && cout << "Block - " << i << " Start - " << start << " End - " << end << endl;
    int j = prime_i; //local multiple

    while(1+2*prime_i <= sqrt(limit*2)) {
      if((1+2*prime_i)*(1+2*prime_i) > end) {
        j = end;
        if(debug_msg_first_index) cout << "Block - " << i << " IGNORED J" << endl;
      } else {
        if(start == 0 ) {
          j = start + prime_i;
        } else {
          j = start + start % prime_i;
        }
        if(debug_msg_first_index) cout << "Block - " << i << " -  First_J - " << j << endl;
      }

      if(debug_msg_primes_used) cout << "Block - " << i << " Current j = " << j << endl;

      for (int k = j + (j * (2 * j + 1)); k <= end; k += 2 * j + 1) {
        if(debug_msg_marking) cout << "Block - " << i << " k - " << k << endl;
        is_prime[k] = false;
      }

      #pragma omp barrier
      if(i == 0) { //choose new prime tester
        prime_i++;
        for(; prime_i  < limit; prime_i++) {
          if(is_prime[prime_i]) {
            j = prime_i;
            if(debug_msg_primes_used) cout << "Block - " << i << " NEW prime_i - " << prime_i << endl << endl;
            break;
          }
        }
      }
      #pragma omp barrier
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

  omp_set_num_threads(omp_get_max_threads()) ;
  n_threads = omp_get_max_threads();

  int limit_t = floor(limit/2.0);

  vector<bool> is_prime(limit_t, true);

  sieveBlockwise(limit_t, is_prime);

  //Print
  cout << "Primes up to " << limit << ":" << endl;
  for(int i = 0; i < is_prime.size(); i++) {
    is_prime[i] && cout << 2*i+1 << " ";
    //cout << is_prime[i] << endl;
  }

  return EXIT_SUCCESS;
}