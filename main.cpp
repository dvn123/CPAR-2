#include <iostream>
#include <limits>
#include <sstream>
#include <cmath>
#include <vector>

using namespace std;

int sieve(int limit, vector<bool>&is_prime) {
  for(int i = 2; i <= sqrt(limit); i++) {
    if(is_prime[i]) {
      for (int j = i * i; j <= limit; j += i) {
        is_prime[j] = false;
      }
    }
  }
}

int main(int argc, char *argv[]) {
  int limit = 10000;
  if (argc == 2) {
    stringstream ss(argv[--argc]);
    ss >> limit;

    if (limit < 1 or ss.fail()) {
      cerr << "USAGE:\n  sieve LIMIT\n\nwhere LIMIT in the range [1, " << numeric_limits<int>::max() << "]" << endl;
      return 2;
    }
  }

  vector<bool> is_prime(limit, true);

  sieve(limit, is_prime);

  //Print
  //cout << "Primes up to " << limit << ":" << endl;
  for(int i = 2; i < is_prime.size(); i++) {
    (is_prime[i] && cout << i << " ");
  }
}