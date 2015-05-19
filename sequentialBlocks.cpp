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

#define REAL_NUMBER(n) \
(1+2*n)

#define ARRAY_INDEX(n) \
((n-1)/2)

using namespace std;

bool debug_msg_primes_used = false;
bool debug_msg_marking = false;
bool debug_msg_block_assignment = false;
bool debug_msg_first_index = false;

int n_threads;

void sieveBlockwise(int limit, vector<bool> &is_prime) {
  for(long thread_id = 0; thread_id <= n_threads - 1; thread_id++) {
    long prime_i = 1;

    long start = BLOCK_LOW(thread_id, n_threads, limit);
    long end = BLOCK_HIGH(thread_id,n_threads,limit);

    if(debug_msg_block_assignment) cout << endl << "Block - " << thread_id << " Start - " << start << ", End - " << end << "; Real Start - " << REAL_NUMBER(start) << ", End - " << REAL_NUMBER(end) << endl;
    long j; //local multiple

    while(REAL_NUMBER(prime_i)*REAL_NUMBER(prime_i) <= limit*2) {
      if(prime_i > start) {
        j = end;
        if(debug_msg_first_index) cout << "Block - " << thread_id << " Ignored J" << endl;
      } else if(REAL_NUMBER(prime_i)*REAL_NUMBER(prime_i) > REAL_NUMBER(end)) {
        j = end;
        if(debug_msg_first_index) cout << "Block - " << thread_id << " Ignored J" << endl;
      } else {
        if(REAL_NUMBER(start) % REAL_NUMBER(prime_i) == 0) {
          j = start;
        } else {
          j = ARRAY_INDEX(REAL_NUMBER(start) - (REAL_NUMBER(start) % REAL_NUMBER(prime_i)) + REAL_NUMBER(prime_i));
        }
        if(debug_msg_first_index) cout << "Block - " << thread_id << " First_J - " << j << "; Real - " << REAL_NUMBER(j) << endl;
      }

      if(debug_msg_primes_used) cout << "Block - " << thread_id << " Current j = " << j << "; Real - " << REAL_NUMBER(j) << endl;

      for (long k = j; k <= end; k += REAL_NUMBER(prime_i)) {
        if(debug_msg_marking) cout << "Block - " << thread_id << " k - " << k << "; Real - " << REAL_NUMBER(k) << endl;
        is_prime[k] = false;
      }

      //choose new prime
      if(debug_msg_primes_used) cout << endl;
      prime_i++;
      for(; prime_i  < limit; prime_i++) {
        if(is_prime[prime_i]) {
          j = prime_i;
          if(debug_msg_primes_used) cout << "Block - " << thread_id << " prime_i - " << prime_i << "; Real - " << REAL_NUMBER(prime_i) << endl;
          break;
        }
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
      return EXIT_FAILURE;
    }
  } else {
    cerr << "USAGE:\n  sieve LIMIT\n\nwhere LIMIT in the range [1, " << numeric_limits<long>::max() << "]" << endl;
    return EXIT_FAILURE;
  }

  omp_set_num_threads(omp_get_max_threads()) ;
  n_threads = omp_get_max_threads();

  long limit_t = floor(limit/2.0);

  vector<bool> is_prime(limit_t, true);

  sieveBlockwise(limit_t, is_prime);

  //Print
  cout << "Primes up to " << limit << ":" << endl << "2 ";
  for(long i = 1; i < is_prime.size(); i++) {
    is_prime[i] && cout << REAL_NUMBER(i) << " ";
    //cout << is_prime[i] << endl;
  }

  return EXIT_SUCCESS;
}