#include <iostream>
#include <limits>
#include <sstream>
#include <cmath>
#include <vector>
#include <omp.h>
#include <sys/time.h>


#define BLOCK_LOW(id, p, n) \
((id)*(n)/(p))

#define BLOCK_HIGH(id, p, n) \
(BLOCK_LOW((id)+1,p,n)-1)

#define REAL_NUMBER(n) \
(2*n+1)

#define ARRAY_INDEX(n) \
((n-1)/2)

using namespace std;

bool debug_msg_primes_used = false;
bool debug_msg_marking = false;
bool debug_msg_block_assignment = false;
bool debug_msg_first_index = false;
bool debug_msg_primes_found = false;

bool debug_msg_primes_seed = true;


vector<int> blocks_allowed_to_print;

void print(string s, int block) {
    for (int i = 0; i < blocks_allowed_to_print.size(); i++) {
        cout << s;
    }
}

int main(int argc, char *argv[]) {
    blocks_allowed_to_print.push_back(0);
    //blocks_allowed_to_print.push_back(1);
    //blocks_allowed_to_print.push_back(2);
    //blocks_allowed_to_print.push_back(3);


    int n_threads = omp_get_max_threads();

    long limit = -1;
    if (argc == 3) {
        stringstream ss(argv[argc - 2]);
        ss >> limit;

        stringstream ss2(argv[argc - 1]);
        ss2 >> n_threads;

        if (n_threads > omp_get_max_threads())
            n_threads = omp_get_max_threads();

        if (limit < 1 || ss.fail()) {
            cerr << "USAGE:\n  sieve LIMIT n_threads\nwhere LIMIT in the range [1, " << numeric_limits<long>::max() <<
            "] " << endl;
            return -1;
        }
    } else {
        cerr << "USAGE:\n  sieve LIMIT n_threads\n\nwhere LIMIT in the range [1, " << numeric_limits<long>::max() <<
        "]" << endl;
        return -1;
    }

    if ((2 + (limit - 1 / n_threads)) < (int) sqrt((double) limit)) {
        if (n_threads == 0) printf("Too many processes.\n");
        return -1;
    }
    long global_count = 0;
    omp_set_num_threads(n_threads);
    long limit_t = floor(limit / 2.0);

    clock_t begin_time = clock();

    struct timeval t_start;
    gettimeofday(&t_start, NULL);
    int thread_id;

    #pragma omp parallel private(thread_id) num_threads(n_threads)
    {
        thread_id = omp_get_thread_num();

        long start = BLOCK_LOW(thread_id, n_threads, limit_t);
        long end = BLOCK_HIGH(thread_id, n_threads, limit_t);
        long end_first = BLOCK_HIGH(0, n_threads, limit_t);
        if (debug_msg_block_assignment)
            print(string("\nBlock - ") + to_string(thread_id) + " Start - " + to_string(start) + ", End - " +
                  to_string(end) + "; Real Start - " + to_string(REAL_NUMBER(start)) + ", End - " +
                  to_string(REAL_NUMBER(end)) + "\n", thread_id);


        vector<bool> is_prime((end - start) + 1, true);
        vector<bool> is_prime_first((end - start)*4 + 1, true);
        vector<long> primes_to_sieve;

        is_prime_first[0] = false;
        for (long i = 1; REAL_NUMBER(i) * REAL_NUMBER(i) < end_first * 2; i++) {
            if (is_prime_first[i]) {
                if (debug_msg_primes_seed) print(string("\n[Seed List] Block - ") + to_string(0) + " i - " + to_string(i) + "; Real - " + to_string(REAL_NUMBER(i)) + "\n", thread_id);
                for (long h = i + (i * (REAL_NUMBER(i))); h * h < end_first; h += REAL_NUMBER(i)) {
                    if (debug_msg_primes_seed) print(string("[Seed List] Block - ") + to_string(0) + " k - " + to_string(h) + "; Real - " + to_string(REAL_NUMBER(h)) + "; Vector - " + to_string(h) + "\n", thread_id);
                    is_prime_first[h] = false;
                }
            }
        }
        if (debug_msg_primes_seed) print(string("Primes list - "), thread_id);

        for (long i = 0; i < sqrt(is_prime_first.size()); i++) {
            if (is_prime_first[i]) {
                if (debug_msg_primes_seed) print(string(to_string(REAL_NUMBER(i))) + " ", thread_id);
                primes_to_sieve.push_back(i);
            }
        }
        if (debug_msg_primes_seed) print(string("\n\n"), thread_id);


        long prime_i = 1;
        long prime_index = 0;

        if (debug_msg_block_assignment)
            print(string("\nBlock - ") + to_string(thread_id) + " Start - " + to_string(start) + ", End - " +
                  to_string(end) + "; Real Start - " + to_string(REAL_NUMBER(start)) + ", End - " +
                  to_string(REAL_NUMBER(end)) + "\n", thread_id);
        long j; //local multiple


        while (REAL_NUMBER(prime_i) * REAL_NUMBER(prime_i) < limit_t * 2) {
            if (REAL_NUMBER(prime_i) * REAL_NUMBER(prime_i) > REAL_NUMBER(end)) {
                j = end + 1;
                if (debug_msg_first_index)
                    print(string("\nBlock - ") + to_string(thread_id) + " Ignored_J \n", thread_id);
            } else {
                if (REAL_NUMBER(start) % REAL_NUMBER(prime_i) == 0) {
                    j = start;
                } else {
                    if ((REAL_NUMBER(start) - (REAL_NUMBER(start) % REAL_NUMBER(prime_i)) + REAL_NUMBER(prime_i)) % 2 ==
                        0) {
                        j = ARRAY_INDEX(REAL_NUMBER(start) - (REAL_NUMBER(start) % REAL_NUMBER(prime_i)) +
                                        2 * REAL_NUMBER(prime_i));
                    } else
                        j = ARRAY_INDEX(REAL_NUMBER(start) - (REAL_NUMBER(start) % REAL_NUMBER(prime_i)) +
                                        REAL_NUMBER(prime_i));
                }
                if (debug_msg_first_index)
                    print(string("Block - ") + to_string(thread_id) + " First_J - " + to_string(j) + "; Real - " +
                          to_string(REAL_NUMBER(j)) + "\n", thread_id);
            }
            if (REAL_NUMBER(j) < REAL_NUMBER(prime_i) * REAL_NUMBER(prime_i)) {
                j = ARRAY_INDEX(REAL_NUMBER(prime_i) * REAL_NUMBER(prime_i));
            }


            for (long k = j; k <= end; k += REAL_NUMBER(prime_i)) {
                if (debug_msg_marking)
                    print(string("Block - ") + to_string(thread_id) + " k - " + to_string(k) + "; Real - " +
                          to_string(REAL_NUMBER(k)) + "; Vector - " + to_string(k - start) + "\n", thread_id);
                if (prime_i != k)
                    is_prime[k - start] = false;
            }
            prime_i = primes_to_sieve[++prime_index];
            if (debug_msg_primes_used)
                print(string("\nBlock - ") + to_string(thread_id) + " prime_i - " + to_string(prime_i) + "; Real - " +
                      to_string(REAL_NUMBER(prime_i)) + "\n", thread_id);
        }

        /*long count = 0;
        for (long i = 0; i < is_prime.size(); i++) {
            //print(to_string(is_prime[i]) + " ", thread_id);
            if (is_prime[i] == true) {
                //if (debug_msg_primes_found) print(to_string(REAL_NUMBER(i) + (REAL_NUMBER(start) - 1)) + ", ", thread_id);
                count++;
            }
        }
        //cout << endl << endl << "Block - " << thread_id << " Count = " << count << endl;

        #pragma omp atomic
        global_count += count;*/
    }

    clock_t end_time = clock();
    struct timeval t_end;
    gettimeofday(&t_end, NULL);

    //Print
    cout << "Counted primes up to " << limit << ": " << global_count << endl;
    cout << "Time: " <<
    (t_end.tv_sec + (t_end.tv_usec * pow(10, -6))) - (t_start.tv_sec + (t_start.tv_usec / pow(10, 6))) << endl;
    cout << "CPU Time: " << (end_time - begin_time) / double(CLOCKS_PER_SEC) << endl;


    return 0;
}