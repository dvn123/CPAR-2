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
vector<int> blocks_allowed_to_print;

void print(string  s, int block) {
  for(int i = 0; i < blocks_allowed_to_print.size(); i++) {
    if(block == blocks_allowed_to_print[i])
      cout << s;
  }
}

int n_threads;

void sieveBlockwise(int limit, vector<bool> &is_prime) {
  int prime_i = 1;

  int thread_n;

  #pragma omp parallel for default(shared) private(thread_n)
  for(thread_n = 0; thread_n <= n_threads - 1; ++thread_n) {
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
          int i_temp = 1;
          while((REAL_NUMBER(start) - (REAL_NUMBER(start) % REAL_NUMBER(prime_i)) + i_temp*REAL_NUMBER(prime_i)) % 2 == 0) {
            i_temp++;
          }
          j = ARRAY_INDEX(REAL_NUMBER(start) - (REAL_NUMBER(start) % REAL_NUMBER(prime_i)) + i_temp*REAL_NUMBER(prime_i));
        }
        if(debug_msg_first_index) cout << "Block - " << thread_n << " First_J - " << j << "; Real - " << REAL_NUMBER(j) << endl;
      }

      if(debug_msg_primes_used) cout << "Block - " << thread_n << " Current j = " << j << "; Real - " << REAL_NUMBER(j) << endl;

      for (int k = j; k <= end; k += REAL_NUMBER(prime_i)) {
        if(debug_msg_marking) print(string("Block - ") + to_string(thread_n) + " k - " + to_string(k) + "; Real - " + to_string(REAL_NUMBER(k)) + "; Vector - " + to_string(k-start) + "\n", thread_n);
        is_prime[k-start] = false;
      }


      //choose new prime
      if(debug_msg_primes_used) cout << endl;
      #pragma omp barrier
      if(thread_n == 0) {
        for(++prime_i; prime_i  < limit; prime_i++) {
          if(is_prime[prime_i]) {
            j = prime_i;
            if(debug_msg_primes_used) cout << "Block - " << thread_n << " prime_i - " << prime_i << "; Real - " << REAL_NUMBER(prime_i) << endl;
            break;
          }
        }
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