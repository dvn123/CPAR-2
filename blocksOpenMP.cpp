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

#define REAL_NUMBER(n) \
(1+2*n)

#define ARRAY_INDEX(n) \
((n-1)/2)

using namespace std;

bool debug_msg_primes_used = true;
bool debug_msg_marking = true;
bool debug_msg_block_assignment = true;
bool debug_msg_first_index = true;

int n_threads;

void sieveBlockwise(int limit, vector<bool> &is_prime) {
  int prime_i = 1;

  #pragma omg parallel for default(shared) private(i)
  for(int thread_n = 0; thread_n <= n_threads - 1; ++thread_n) {
    int start = BLOCK_LOW(thread_n, n_threads, limit);
    int end = BLOCK_HIGH(thread_n,n_threads,limit);

    if(debug_msg_block_assignment) cout << endl << "Block - " << thread_n << " Start - " << start << ", End - " << end << "; Real Start - " << REAL_NUMBER(start) << ", End - " << REAL_NUMBER(end) << endl;
    int j; //local multiple

    while(REAL_NUMBER(prime_i)*REAL_NUMBER(prime_i) <= limit*2) {
      if(prime_i > end) {
        j = end;
        if(debug_msg_first_index) cout << "Block - " << thread_n << " Ignored J" << endl;
      } else if(REAL_NUMBER(prime_i)*REAL_NUMBER(prime_i) > REAL_NUMBER(end)) {
        j = end;
        if(debug_msg_first_index) cout << "Block - " << thread_n << " Ignored J" << endl;
      } else {
        if(REAL_NUMBER(start) % REAL_NUMBER(prime_i) == 0) {
          j = start;
        } else {
          j = ARRAY_INDEX(REAL_NUMBER(start) - (REAL_NUMBER(start) % REAL_NUMBER(prime_i)) + REAL_NUMBER(prime_i));
        }
        if(debug_msg_first_index) cout << "Block - " << thread_n << " First_J - " << j << "; Real - " << REAL_NUMBER(j) << endl;
      }

      if(debug_msg_primes_used) cout << "Block - " << thread_n << " Current j = " << j << "; Real - " << REAL_NUMBER(j) << endl;

      for (int k = j; k <= end; k += REAL_NUMBER(prime_i)) {
        if(debug_msg_marking) cout << "Block - " << thread_n << " k - " << k << "; Real - " << REAL_NUMBER(k) << endl;
        is_prime[k] = false;
      }
      #pragma omp barrier


      //choose new prime
      if(debug_msg_primes_used) cout << endl;
      #pragma omp master
      {
        for(++prime_i; prime_i  < limit; prime_i++) {
          if(is_prime[prime_i]) {
            j = prime_i;
            if(debug_msg_primes_used) cout << "Block - " << thread_n << " prime_i - " << prime_i << "; Real - " << REAL_NUMBER(prime_i) << endl;
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

  n_threads = omp_get_max_threads();

  omp_set_num_threads(n_threads) ;

  int limit_t = floor(limit/2.0);

  vector<bool> is_prime(limit_t, true);

  sieveBlockwise(limit_t, is_prime);

  //Print
  cout << "Primes up to " << limit << ":" << endl << "2 ";
  for(int i = 1; i < is_prime.size(); i++) {
    is_prime[i] && cout << REAL_NUMBER(i) << " ";
  }

  return EXIT_SUCCESS;
}