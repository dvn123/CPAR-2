#include <iostream>
#include <limits>
#include <sstream>
#include <cmath>
#include <vector>
#include <mpi.h>

#define EXIT_SUCCESS 0
#define EXIT_FAILURE -1

#define BLOCK_LOW(id,p,n) \
((id)*(n)/(p))

#define BLOCK_HIGH(id,p,n) \
(BLOCK_LOW((id)+1,p,n)-1)

#define REAL_NUMBER(n) \
(2*n+1)

#define ARRAY_INDEX(n) \
((n-1)/2)

#define REAL_LIMIT(n) \
(n*2)

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
vector<int> blocks_allowed_to_print;

int thread_id;

void print(string  s, bool type) {
  if(!type)
    return;
  for(int i = 0; i < blocks_allowed_to_print.size(); i++) {
    if(thread_id == blocks_allowed_to_print[i])
      cout << s;
  }
}

int checkArgs(int argc, char *argv[], int &limit, int n_threads) {
  blocks_allowed_to_print.push_back(0);

  if (argc == 2) {
    stringstream ss(argv[--argc]);
    ss >> limit;
    limit = floor(limit/2.0);

    if (limit < 1 || ss.fail()) {
      cerr << "USAGE:\n  sieve LIMIT\nwhere LIMIT in the range [1, " << numeric_limits<int>::max() << "]" << endl;
      return EXIT_FAILURE;
    }
  } else {
    cerr << "USAGE:\n  sieve LIMIT\n\nwhere LIMIT in the range [1, " << numeric_limits<int>::max() << "]" << endl;
    return EXIT_FAILURE;
  }

  if ((2 + (REAL_LIMIT(limit) - 1 / n_threads)) < (int) sqrt((double) REAL_LIMIT(limit))) {
    if (n_threads == 0) printf("Too many processes.\n");
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

int getFirstMultipleInRange(int prime, int start, int end) {
  if(REAL_NUMBER(prime)*REAL_NUMBER(prime) > REAL_NUMBER(end)) {
    print(string("\nBlock - ") + to_string(thread_id) + " Ignored - " + to_string(prime) + "; Real -" + to_string(REAL_NUMBER(prime)) + "\n", debug_msg_first_index);
    return end + 1;
  } else {
    if(REAL_NUMBER(start) % REAL_NUMBER(prime) == 0) {
      return start;
    } else {
      int i_temp = 1;
      while((REAL_NUMBER(start) - (REAL_NUMBER(start) % REAL_NUMBER(prime)) + i_temp*REAL_NUMBER(prime)) % 2 == 0) {
        i_temp++;
      }
      return ARRAY_INDEX(REAL_NUMBER(start) - (REAL_NUMBER(start) % REAL_NUMBER(prime)) + i_temp*REAL_NUMBER(prime));
    }
  }
}

int chooseNewPrime(int prime, int limit, vector<bool> &is_prime) {
  for(prime = ++prime; prime  < limit; prime++) {
    //cout << is_prime[prime] << endl;
    if (is_prime[prime]) {
      print(string("\nBlock - ") + to_string(thread_id) + " original_prime - " + to_string(prime) + "; Real - " + to_string(REAL_NUMBER(prime)) + "\n", debug_msg_primes_used);
      return prime;
    }
  }
  MPI_Finalize();
  return EXIT_FAILURE;
}

void sieveBlock(int start, int end, int limit, vector<bool> &is_prime, int end_first, vector<bool> &is_prime_first) {
  int original_prime = 1;

  int local_multiple;

  while (REAL_NUMBER(original_prime) * REAL_NUMBER(original_prime) <= REAL_LIMIT(limit)) {
    local_multiple = getFirstMultipleInRange(original_prime, start, end);
    print(string("Block - ") + to_string(thread_id) + " First Multiple - " + to_string(local_multiple) + "; Real - " +
                                                                                                         to_string(REAL_NUMBER(local_multiple)) + "\n", debug_msg_first_index);

    //Mark the multiples in the block
    for (int k = local_multiple; k <= end; k += REAL_NUMBER(original_prime)) {
      print(string("Block - ") + to_string(thread_id) + " k - " + to_string(k) + "; Real - " + to_string(REAL_NUMBER(k)) + "; Vector - " + to_string(k - start) + "\n", debug_msg_marking);
      if(original_prime != k)
          is_prime[k-start] = false;
    }

    //Every thread has its own list - improvement 2
    if(thread_id != 0) {
      for (int k = original_prime; k <= end_first; k += REAL_NUMBER(original_prime)) {
        print(string("Block - ") + to_string(thread_id) + " fake k - " + to_string(k) + "; Real - " + to_string(REAL_NUMBER(k)) + "; Vector - " + to_string(k) + "\n", debug_msg_marking);
        if (original_prime != k)
          is_prime_first[k] = false;
      }
    }

    if(thread_id != 0)
      original_prime = chooseNewPrime(original_prime, limit, is_prime_first);
    else
      original_prime = chooseNewPrime(original_prime, limit, is_prime);
  }
}

void sieve(int n_threads, int limit, int &count) {
  double runtime = -MPI_Wtime();

  int local_count;
  int start = BLOCK_LOW(thread_id, n_threads, limit);
  int end = BLOCK_HIGH(thread_id, n_threads, limit);
  vector<bool> is_prime(end-start+1, true);

  int end_first = BLOCK_HIGH(0, n_threads, limit);
  vector<bool> is_prime_first;

  if(thread_id != 0) {
    for(int i = 0; i < end_first; i++)
      is_prime_first.push_back(true);
  }


  local_count = 0;
  print(string("\nBlock - ") + to_string(thread_id) + " Start - " + to_string(start) + ", End - " + to_string(end) + "; Real Start - " + to_string(REAL_NUMBER(start)) + ", End - " + to_string(REAL_NUMBER(end)) + "\n", debug_msg_block_assignment);

  sieveBlock(start, end, limit, is_prime, end_first, is_prime_first);

  for (int i = 0; i < is_prime.size(); i++) {
    if (is_prime[i] == true) {
      print(to_string(REAL_NUMBER(i) + (REAL_NUMBER(start)-1)) + ", ", debug_msg_primes_found);
      local_count++;
    }
  }

  MPI_Reduce(&local_count, &count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  runtime += MPI_Wtime();

  //Print
  if (thread_id == 0) {
    cout << "Counted primes up to " << REAL_LIMIT(limit) << ": " << count << endl << "Time: " << runtime << endl;
  }
}

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  int limit = -1, n_threads = 1;
  int count = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &thread_id);
  MPI_Comm_size(MPI_COMM_WORLD, &n_threads);

  if(checkArgs(argc, argv, limit, n_threads) != EXIT_SUCCESS) {
    MPI_Finalize();
    return EXIT_FAILURE;
  }

  sieve(n_threads, limit, count);

  MPI_Finalize();
  return 0;
}