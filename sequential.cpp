#include <iostream>
#include <limits>
#include <sstream>
#include <cmath>
#include <vector>
#include <ctime>
#include <sys/time.h>

#define REAL_NUMBER(n) \
(1+2*n)

#define ARRAY_INDEX(n) \
((n-1)/2)

using namespace std;

bool debug_msg_i = false;
bool debug_msg_j = false;

void sieve(long limit, vector<bool> &is_prime) {
  for(long i = 1; REAL_NUMBER(i)*REAL_NUMBER(i) <= limit*2; i++) {
    debug_msg_i && cout << "i - " << i << endl;
    debug_msg_i && cout << REAL_NUMBER(i) << " < " << sqrt(limit*2) << endl;
    if(is_prime[i]) {
      debug_msg_i && cout << "Chosen" << endl;
      for (long j = i + (i*(REAL_NUMBER(i))); j <= limit; j += REAL_NUMBER(i)) {
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
      return EXIT_FAILURE;
    }
  } else {
    cerr << "USAGE:\n  sieve LIMIT\n\nwhere LIMIT in the range [1, " << numeric_limits<long>::max() << "]" << endl;
    return EXIT_FAILURE;
  }

  long limit_t = floor(limit/2.0);

  vector<bool> is_prime(limit_t, true);

  timespec start, end;

  struct timeval t_start;
  gettimeofday(&t_start,NULL);

  (CLOCK_REALTIME, &start);
  clock_t begin_cpu = clock();


  sieve(limit_t, is_prime);

  clock_t end_cpu = clock();
  clock_gettime(CLOCK_REALTIME, &end);

  struct timeval t_end;
  gettimeofday(&t_end,NULL);


  long count = 0;

  for(long i = 0; i < is_prime.size(); i++) {
    is_prime[i] && count++;
  }
  //Print
  cout << "Number of primes up to " << limit << ": " << count << endl;
  cout << "Time: " << (t_end.tv_sec + (t_end.tv_usec*pow(10,-6)))-(t_start.tv_sec+(t_start.tv_usec/pow(10,6))) << endl;
  cout << "CPU Time: " << (end_cpu-begin_cpu)/double(CLOCKS_PER_SEC) << endl;

  return EXIT_SUCCESS;
}