#include <iostream>
#include <limits>
#include <sstream>
#include <cmath>
#include <vector>
#include <omp.h>
#include <mpi.h>

#define BLOCK_LOW(id,p,n) \
((id)*(n)/(p))

#define BLOCK_HIGH(id,p,n) \
(BLOCK_LOW((id)+1,p,n)-1)

#define REAL_NUMBER(n) \
(2*n+1)

#define ARRAY_INDEX(n) \
((n-1)/2)

//mpi doesn't support c++11...
template < typename T > std::string to_string( const T& n )
{
  std::ostringstream stm ;
  stm << n ;
  return stm.str() ;
}

using namespace std;

bool debug_msg_primes_used = false;
bool debug_msg_marking = false;
bool debug_msg_block_assignment = false;
bool debug_msg_first_index = false;
bool debug_msg_primes_found = false;

bool debug_msg_primes_seed = false;

vector<int> blocks_allowed_to_print;

int n_threads, thread_id;

void print(string  s, int block) {
  for(int i = 0; i < blocks_allowed_to_print.size(); i++) {
    if(block == blocks_allowed_to_print[i])
      cout << s;
  }
}

void sieveBlockwise(long limit, vector<bool> &is_prime, long start, long end, vector<bool> &is_prime_first, long end_first, vector<long> primes_to_sieve) {
  long prime_i = 1;
    long prime_index = 0;

  if(debug_msg_block_assignment) print(string("\nBlock - ") + to_string(thread_id) + " Start - " + to_string(start) + ", End - " + to_string(end) + "; Real Start - " + to_string(REAL_NUMBER(start)) + ", End - " + to_string(REAL_NUMBER(end)) + "\n", thread_id);
  long j; //local multiple



  while(REAL_NUMBER(prime_i)*REAL_NUMBER(prime_i) <= limit*2) {

    if(REAL_NUMBER(prime_i)*REAL_NUMBER(prime_i) > REAL_NUMBER(end)) {
      j = end+1;
      if(debug_msg_first_index)  print(string("\nBlock - ") + to_string(thread_id) + " Ignored_J \n", thread_id);
    } else {
      if(REAL_NUMBER(start) % REAL_NUMBER(prime_i) == 0) {
        j = start;
      } else {
        long i_temp = 1;
        while((REAL_NUMBER(start) - (REAL_NUMBER(start) % REAL_NUMBER(prime_i)) + i_temp*REAL_NUMBER(prime_i)) % 2 == 0) {
          i_temp++;
        }
        j = ARRAY_INDEX(REAL_NUMBER(start) - (REAL_NUMBER(start) % REAL_NUMBER(prime_i)) + i_temp*REAL_NUMBER(prime_i));
      }
      if(debug_msg_first_index) print(string("Block - ") + to_string(thread_id) + " First_J - " + to_string(j) + "; Real - " + to_string(REAL_NUMBER(j)) + "\n", thread_id);
    }
    if(REAL_NUMBER(j) < REAL_NUMBER(prime_i)*REAL_NUMBER(prime_i)) {
      j = ARRAY_INDEX(REAL_NUMBER(prime_i)*REAL_NUMBER(prime_i));
    }


    for (long k = j; k <= end; k += REAL_NUMBER(prime_i)) {
      if(debug_msg_marking) print(string("Block - ") + to_string(thread_id) + " k - " + to_string(k) + "; Real - " + to_string(REAL_NUMBER(k)) + "; Vector - " + to_string(k-start) + "\n", thread_id);
      if(prime_i != k)
        is_prime[k-start] = false;
    }

    //choose new prime
    prime_i = primes_to_sieve[++prime_index];

    //MPI_Bcast(&prime_i,  1, MPI_INT, 0, MPI_COMM_WORLD);
    if (debug_msg_primes_used) print(string("\nBlock - ") + to_string(thread_id) + " prime_i - " + to_string(prime_i) + "; Real - " + to_string(REAL_NUMBER(prime_i)) + "\n", thread_id);
  }
}


int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  blocks_allowed_to_print.push_back(3);
  //blocks_allowed_to_print.push_back(1);

  MPI_Comm_rank(MPI_COMM_WORLD, &thread_id);
  MPI_Comm_size(MPI_COMM_WORLD, &n_threads);

  double runtime = -MPI_Wtime();

  long limit = -1;
  if (argc == 2) {
    stringstream ss(argv[--argc]);
    ss >> limit;

    if (limit < 1 || ss.fail()) {
      cerr << "USAGE:\n  sieve LIMIT\nwhere LIMIT in the range [1, " << numeric_limits<int>::max() << "]" << endl;
      MPI_Finalize();
      return -1;
    }
  } else {
    cerr << "USAGE:\n  sieve LIMIT\n\nwhere LIMIT in the range [1, " << numeric_limits<int>::max() << "]" << endl;
    MPI_Finalize();
    return -1;
  }

  if ((2 + (limit - 1 / n_threads)) < (int) sqrt((double) limit)) {
    if (n_threads == 0) printf("Too many processes.\n");
    MPI_Finalize();
    return -1;
  }

  long limit_t = floor(limit/2.0);

  long start = BLOCK_LOW(thread_id, n_threads, limit_t);
  long end = BLOCK_HIGH(thread_id,n_threads,limit_t);

  long end_first = BLOCK_HIGH(0, n_threads, limit_t);

  vector<bool> is_prime((end - start) + 1, true);
  vector<bool> is_prime_first((end - start) + 1, true);
  vector<long> primes_to_sieve;


  is_prime_first[0]= false;
  for (long i = 1; REAL_NUMBER(i)*REAL_NUMBER(i) <= limit_t*2; i++) {
    if(is_prime_first[i]) {
      if (debug_msg_primes_seed) print(string("\n[Seed List] Block - ") + to_string(thread_id) + " i - " + to_string(i) + "; Real - " + to_string(REAL_NUMBER(i)) + "\n", thread_id);
      for (long h = i + (i * (REAL_NUMBER(i))); h <= sqrt(limit_t); h += REAL_NUMBER(i)) {
        if(debug_msg_primes_seed) print(string("[Seed List] Block - ") + to_string(thread_id) + " k - " + to_string(h) + "; Real - " + to_string(REAL_NUMBER(h)) + "; Vector - " + to_string(h-start) + "\n", thread_id);
        is_prime_first[h] = false;
      }
    }
  }
  if (debug_msg_primes_seed) print(string("Primes list - "), thread_id);
  for (long i = 0; i <= ARRAY_INDEX(sqrt(limit_t*2)); i++) {
    if(is_prime_first[i]) {
      primes_to_sieve.push_back(i);
      if (debug_msg_primes_seed) print(string(to_string(REAL_NUMBER(i))) + " ", thread_id);
    }
  }
  if (debug_msg_primes_seed) print(string("\n\n"), thread_id);

  is_prime_first.clear();


  sieveBlockwise(limit_t, is_prime, start, end, is_prime_first, end_first, primes_to_sieve); // sieveBlockwise(long limit, vector<bool> &is_prime, long start, long end, vector<bool> &is_prime_first, long end_first) {


  /*long count = 0;
  for (long i = 0; i < is_prime.size(); i++) {
    //print(to_string(is_prime[i]) + " ", thread_id);
    if (is_prime[i] == true) {
      if(debug_msg_primes_found) print(to_string(REAL_NUMBER(i) + (REAL_NUMBER(start)-1)) + ", ", thread_id);
      count++;
    }
  }*/

  //print(to_string(is_prime[123]), thread_id);

  long global_count;

  //MPI_Reduce(&count, &global_count, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
  runtime += MPI_Wtime();

  //Print
  if(thread_id == 0) {
    cout << "\nCounted primes up to " << limit << ": "  << global_count << endl << "Time: " << runtime << endl;
  }

  MPI_Finalize();
  return 0;
}