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
bool debug_msg_block_assignment = true;
bool debug_msg_first_index = false;
bool debug_msg_primes_found = false;
vector<int> blocks_allowed_to_print;

int n_threads, thread_id;

void print(string  s, int block) {
    for(int i = 0; i < blocks_allowed_to_print.size(); i++) {
        if(block == blocks_allowed_to_print[i])
            cout << s;
    }
}

int main(int argc, char *argv[]) {

    blocks_allowed_to_print.push_back(0);
    //blocks_allowed_to_print.push_back(1);
    n_threads = omp_get_max_threads();

    int limit = -1;
    if (argc == 3) {
        stringstream ss(argv[argc-2]);
        ss >> limit;

        stringstream ss2(argv[argc-1]);
        ss2 >> n_threads;

        if(n_threads > omp_get_max_threads())
            n_threads = omp_get_max_threads();

        if (limit < 1 || ss.fail()) {
            cerr << "USAGE:\n  sieve LIMIT n_threads\nwhere LIMIT in the range [1, " << numeric_limits<int>::max() << "] " << endl;
            return -1;
        }
    } else {
        cerr << "USAGE:\n  sieve LIMIT n_threads\n\nwhere LIMIT in the range [1, " << numeric_limits<int>::max() << "]" << endl;
        return -1;
    }

    if ((2 + (limit - 1 / n_threads)) < (int) sqrt((double) limit)) {
        if (n_threads == 0) printf("Too many processes.\n");
        return -1;
    }
    int global_count = 0;

    omp_set_num_threads(n_threads);

    clock_t begin_time = clock();



    #pragma omp parallel private(thread_id) num_threads(n_threads)
    {
        thread_id = omp_get_thread_num();
        int limit_t = floor(limit / 2.0);

        int start = BLOCK_LOW(thread_id, n_threads, limit_t);
        int end = BLOCK_HIGH(thread_id, n_threads, limit_t);

        int end_first = BLOCK_HIGH(0, n_threads, limit_t);


        vector<bool> is_prime((end - start) + 1, true);
        vector<bool> is_prime_first((end - start) + 1, true);


        int prime_i = 1;

        if(debug_msg_block_assignment) print(string("\nBlock - ") + to_string(thread_id) + " Start - " + to_string(start) + ", End - " + to_string(end) + "; Real Start - " + to_string(REAL_NUMBER(start)) + ", End - " + to_string(REAL_NUMBER(end)) + "\n", thread_id);
        int j; //local multiple


        while(REAL_NUMBER(prime_i)*REAL_NUMBER(prime_i) <= limit*2) {
            if(REAL_NUMBER(prime_i)*REAL_NUMBER(prime_i) > REAL_NUMBER(end)) {
                j = end+1;
                if(debug_msg_first_index)  print(string("\nBlock - ") + to_string(thread_id) + " Ignored_J \n", thread_id);
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
                if(debug_msg_first_index) print(string("Block - ") + to_string(thread_id) + " First_J - " + to_string(j) + "; Real - " + to_string(REAL_NUMBER(j)) + "\n", thread_id);
            }

            for (int k = j; k <= end; k += REAL_NUMBER(prime_i)) {
                if(debug_msg_marking) print(string("Block - ") + to_string(thread_id) + " k - " + to_string(k) + "; Real - " + to_string(REAL_NUMBER(k)) + "; Vector - " + to_string(k-start) + "\n", thread_id);
                if(prime_i != k)
                    is_prime[k-start] = false;
            }

            for (int h = prime_i; h <= end_first; h += REAL_NUMBER(prime_i)) {
                if(debug_msg_marking) print(string("Block - ") + to_string(thread_id) + " fake k - " + to_string(h) + "; Real - " + to_string(REAL_NUMBER(h)) + "; Vector - " + to_string(h) + "\n", thread_id);
                if(prime_i != h)
                    is_prime_first[h] = false;
            }

            //choose new prime
            if(thread_id == 0) {
                //if (debug_msg_primes_used) cout << endl;
                for (++prime_i; prime_i < limit; prime_i++) {
                    if (is_prime[prime_i]) {
                        j = prime_i;
                        //if (debug_msg_primes_used) print(string("Block - ") + to_string(thread_id) + " prime_i - " + to_string(prime_i) + "; Real - " + to_string(REAL_NUMBER(prime_i)) + "\n\n", thread_id);
                        break;
                    }
                }
            } else {
                for (++prime_i; prime_i < limit; prime_i++) {
                    //print(string("Block - ") + to_string(thread_id) + " prime_i TEST - " + to_string(prime_i) + "; Real - " + to_string(REAL_NUMBER(prime_i)) + "\n", thread_id);
                    if (is_prime_first[prime_i]) {
                        //if (debug_msg_primes_used) print(string("Block - ") + to_string(thread_id) + " prime_i - " + to_string(prime_i) + "; Real - " + to_string(REAL_NUMBER(prime_i)) + "\n\n", thread_id);
                        break;
                    }
                }
            }
            if (debug_msg_primes_used) print(string("\nBlock - ") + to_string(thread_id) + " prime_i - " + to_string(prime_i) + "; Real - " + to_string(REAL_NUMBER(prime_i)) + "\n", thread_id);
        }

        int count = 0;
        for (int i = 0; i < is_prime.size(); i++) {
            //print(to_string(is_prime[i]) + " ", thread_id);
            if (is_prime[i] == true) {
                //if (debug_msg_primes_found) print(to_string(REAL_NUMBER(i) + (REAL_NUMBER(start) - 1)) + ", ", thread_id);
                count++;
            }
        }
       //cout << endl << endl << "Block - " << thread_id << " Count = " << count << endl;

        #pragma omp atomic
        global_count += count;
    }

    clock_t end_time = clock();


    //Print
    if(thread_id == 0) {
        cout << "Counted primes up to " << limit << ": "  << global_count << endl << "Time: " << (end_time-begin_time)/CLOCKS_PER_SEC << endl;
    }

    return 0;
}