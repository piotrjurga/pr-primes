#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <omp.h>

#include <mmintrin.h>
#include <xmmintrin.h>

typedef uint8_t  u8;
typedef uint16_t u16;
typedef uint32_t u32;
typedef uint64_t u64;

typedef int8_t  s8;
typedef int16_t s16;
typedef int32_t s32;
typedef int64_t s64;

#define min(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))

u32 find_primes_seq(u32 begin, u32 end, u32 *primes, u32 *primes_begin) {
    u32 prime_count = 0;

    u32 sqrt_end = sqrt(end);
    if (begin > 2) {
        prime_count = find_primes_seq(2, min(begin, sqrt_end), primes, primes_begin);
        *primes_begin = prime_count;
    } else {
        *primes_begin = 0;
    }

    for (u32 i = begin; i <= end; i++) {
        bool is_prime = true;
        for (u32 j = 0; j < prime_count; j++) {
            u32 p = primes[j];
            if (p*p > i) break;
            if (i % p == 0) {
                is_prime = false;
                break;
            }
        }
        if (is_prime) {
            primes[prime_count++] = i;
        }
    }

    return prime_count;
}

void find_primes_seq_v2(u32 begin, u32 end, u8 *primes) {
    u32 sqrt_end = sqrt(end);
    u32 *base_primes = ((u32 *)primes) + end/4 + 2;
    u32 ignore;
    u32 base_prime_count = find_primes_seq(2, sqrt_end, base_primes, &ignore);

    for (u32 i = begin; i <= end; i++) {
        for (u32 j = 0; j < base_prime_count; j++) {
            u32 p = base_primes[j];
            if (p*p > i) break;
            if (i % p == 0) {
                primes[i] = 1;
                break;
            }
        }
    }
}

#define DEBUG_PARALLEL

u32 find_primes_parallel(u32 begin, u32 end, u32 *primes) {
#ifdef DEBUG_PARALLEL
    auto t1 = omp_get_wtime();
#endif
    u32 prime_count = 0;

    u32 sqrt_end = sqrt(end);
    u32 ignore;
    prime_count = find_primes_seq(2, sqrt_end, primes, &ignore);
#ifdef DEBUG_PARALLEL
    auto t2 = omp_get_wtime();
#endif

#pragma omp parallel shared(prime_count)
    {
#pragma omp for
        for (s32 i = max(begin, sqrt_end+1); i <= end; i++) {
            bool is_prime = true;
            for (u32 j = 0; j < prime_count; j++) {
                u32 p = primes[j];
                if (p*p > i) break;
                if (i % p == 0) {
                    is_prime = false;
                    break;
                }
            }
            if (is_prime) {
#pragma omp critical
                primes[prime_count++] = i;
            }
        }
    }
#ifdef DEBUG_PARALLEL
    auto t3 = omp_get_wtime();
    printf("sequential %fs of all %f (%f%%)\n", t2-t1, t3-t1, 100.0*(t2-t1)/(t3-t1));
#endif

    return prime_count;
}

u32 find_primes_parallel_v3(u32 begin, u32 end, u32 *primes) {
    u32 prime_count = 0;

    u32 sqrt_end = sqrt(end);
    u32 ignore;
    prime_count = find_primes_seq(2, sqrt_end, primes, &ignore);

#pragma omp parallel shared(prime_count)
    {
        u32 local_prime_count = 0;
        u32 *local_primes = (u32 *)malloc(end-begin);
#pragma omp for
        for (s32 i = max(begin, sqrt_end+1); i <= end; i++) {
            bool is_prime = true;
            for (u32 j = 0; j < prime_count; j++) {
                u32 p = primes[j];
                if (p*p > i) break;
                if (i % p == 0) {
                    is_prime = false;
                    break;
                }
            }
            if (is_prime) {
                local_primes[local_prime_count++] = i;
            }
        }
#pragma omp critical
        {
            for (u32 i = 0; i < local_prime_count; i++) {
                primes[prime_count++] = local_primes[i];
            }
        }
        free(local_primes);
    }


    return prime_count;
}

void find_primes_parallel_v2(u32 begin, u32 end, u8 *primes) {
    u32 sqrt_end = sqrt(end);
    u32 *base_primes = ((u32 *)primes) + end/4 + 2;
    u32 ignore;
    u32 base_prime_count = find_primes_seq(2, sqrt_end, base_primes, &ignore);
    end++;

#pragma omp parallel
    {
        u32 work_per_thread = (end-begin) / omp_get_num_threads() + 1;
        u32 thread_begin = begin + work_per_thread*omp_get_thread_num();
        u32 thread_end = thread_begin + work_per_thread;
        if (thread_end > end) thread_end = end;

#if 0
        printf("thread %d starts at %d, ends at %d\n", omp_get_thread_num(), thread_begin, thread_end-1);
#endif

        for (u32 i = thread_begin; i <= thread_end; i++) {
            for (u32 j = 0; j < base_prime_count; j++) {
                u32 p = base_primes[j];
                if (p*p > i) break;
                if (i % p == 0) {
                    primes[i] = 1;
                    break;
                }
            }
        }
    }
}

void find_primes_sieve_seq(u32 begin, u32 end, u8 *primes) {
    u32 sqrt_end = sqrt(end);

    if (begin == 2) {
        for (u32 i = 2; i <= sqrt_end; i++) {
            if (!primes[i]) {
                for (u32 j = i+i; j <= end; j += i) {
                    primes[j] = 1;
                }
            }
        }
    } else {
        find_primes_sieve_seq(2, sqrt_end, primes);
        for (u32 i = 2; i <= sqrt_end; i++) {
            if (!primes[i]) {
                u32 start = begin-1 + i - ((begin-1) % i);
                for (u32 j = start; j <= end; j += i) {
                    primes[j] = 1;
                }
            }
        }
    }
}

void find_primes_sieve_domain(u32 begin, u32 end, u8 *primes) {
    u32 sqrt_end = sqrt(end);

    find_primes_sieve_seq(2, sqrt_end, primes);
    begin = max(sqrt_end+1, begin);

#pragma omp parallel
    {
        u32 work_per_thread = (end-begin) / omp_get_num_threads() + 1;
        u32 thread_begin = begin + work_per_thread*omp_get_thread_num();
        u32 thread_end = thread_begin + work_per_thread;
        if (thread_end > end) thread_end = end;

        for (u32 i = 2; i <= sqrt_end; i++) {
            if (!primes[i]) {
                u32 start = thread_begin-1 + i - ((thread_begin-1) % i);
                for (u32 j = start; j <= thread_end; j += i) {
                    primes[j] = 1;
                }
            }
        }
    }
}

void find_primes_sieve_functional(u32 begin, u32 end, u8 *primes) {
    u32 sqrt_end = sqrt(end);
    u32 *base_primes = (u32 *)&primes[end+1];
    find_primes_sieve_seq(2, sqrt_end, primes);
    u32 base_prime_count = 0;
    for (u32 i = 2; i <= sqrt_end; i++) {
        if (!primes[i]) {
            base_primes[base_prime_count++] = i;
        }
    }

#pragma omp parallel
    {
        u32 thread_begin = omp_get_thread_num();
        u32 thread_end = base_prime_count;

        for (u32 i = thread_begin; i < thread_end; i += omp_get_num_threads()) {
            u32 prime = base_primes[i];
            u32 start = begin-1 + prime - ((begin-1) % prime);
            start = max(start, prime+prime);
            for (u32 j = start; j <= end; j += prime) {
                primes[j] = 1;
            }
        }
    }
}

void find_primes_sieve_seq_v2(u32 begin, u32 end, u8 *primes) {
    u32 sqrt_end = sqrt(end);

    if (begin == 2) {
        for (u32 i = 3; i <= sqrt_end; i+=2) {
            if (!primes[i>>1]) {
                for (u32 j = i+i+i; j <= end; j += i+i) {
                    primes[j>>1] = 1;
                }
            }
        }
    } else {
        find_primes_sieve_seq_v2(2, sqrt_end, primes);
        for (u32 i = 3; i <= sqrt_end; i+=2) {
            if (!primes[i>>1]) {
                u32 start = begin-1 + i - ((begin-1) % i);
                start += i * ((start & 1) == 0);
                for (u32 j = start; j <= end; j += i+i) {
                    primes[j>>1] = 1;
                }
            }
        }
    }
}

#define get_bit(a, i) ((a)[(i)>>3] & (1 << ((i)&7)))
#define set_bit(a, i) ((a)[(i)>>3] |= (1 << ((i)&7)))

void find_primes_sieve_seq_v3(u32 begin, u32 end, u8 *primes) {
    u32 sqrt_end = sqrt(end);

    if (begin == 2) {
        for (u32 i = 2; i <= sqrt_end; i++) {
            if (!get_bit(primes, i)) {
                for (u32 j = i+i; j <= end; j += i) {
                    set_bit(primes, j);
                }
            }
        }
    } else {
        find_primes_sieve_seq_v3(2, sqrt_end, primes);
        for (u32 i = 2; i <= sqrt_end; i++) {
            if (!get_bit(primes, i)) {
                u32 start = begin-1 + i - ((begin-1) % i);
                for (u32 j = start; j <= end; j += i) {
                    set_bit(primes, j);
                }
            }
        }
    }
}

void find_primes_sieve_seq_v4(u32 begin, u32 end, u8 *primes) {
    u32 sqrt_end = sqrt(end);

    if (begin == 2) {
        for (u32 i = 3; i <= sqrt_end; i+=2) {
            if (!get_bit(primes, i>>1)) {
                for (u32 j = i+i+i; j <= end; j += i+i) {
                    set_bit(primes, j>>1);
                }
            }
        }
    } else {
        find_primes_sieve_seq_v4(2, sqrt_end, primes);
        for (u32 i = 3; i <= sqrt_end; i+=2) {
            if (!get_bit(primes, i>>1)) {
                u32 start = begin-1 + i - ((begin-1) % i);
                start += i * ((start & 1) == 0);
                for (u32 j = start; j <= end; j += i+i) {
                    set_bit(primes, j>>1);
                }
            }
        }
    }
}

void find_primes_sieve_seq_v4_no_init(u32 begin, u32 end, u8 *primes) {
    u32 sqrt_end = sqrt(end);

    for (u32 i = 3; i <= sqrt_end; i+=2) {
        if (!get_bit(primes, i>>1)) {
            u32 start = i * ((begin-1)/i + 1);
            start += i * ((start & 1) == 0);
            for (u32 j = start; j <= end; j += i+i) {
                set_bit(primes, j>>1);
                _mm_prefetch((char *)(primes+((j+32*i)>>4)), _MM_HINT_T0);
            }
        }
    }
}

//#define DEBUG_SIEVE_V2

void find_primes_sieve_domain_v2(u32 begin, u32 end, u8 *primes) {
    u32 sqrt_end = sqrt(end);

    find_primes_sieve_seq_v4(2, sqrt_end, primes);
    begin = max(sqrt_end+1, begin);
    begin ^= (begin & 7);
    end++;

#pragma omp parallel
    {
        u32 work_per_thread = (end-begin) / omp_get_num_threads() + 1;
        work_per_thread ^= (work_per_thread & 7);
        work_per_thread += 8;
        u32 thread_begin = begin + work_per_thread*omp_get_thread_num();
        u32 thread_end = thread_begin + work_per_thread;
        if (thread_end > end) thread_end = end;
#ifdef DEBUG_SIEVE_V2
        printf("thread %d starts at %d, ends at %d\n", omp_get_thread_num(), thread_begin, thread_end-1);
        auto start_time = omp_get_wtime();
#endif

        find_primes_sieve_seq_v4_no_init(thread_begin, thread_end-1, primes);

#ifdef DEBUG_SIEVE_V2
        auto end_time = omp_get_wtime();
        printf("thread %d took %f seconds\n", omp_get_thread_num(), end_time-start_time);
#endif
    }
}


#include "test.cpp"

int main() {
    void *memory = malloc(1024*1024*1024);
    memset(memory, 0, 1024*1024*1024);
    u8 *primes = (u8 *)memory;

    u32 ignore;
    u32 count = find_primes_seq(2, 500, (u32 *)primes, &ignore);
    print_list((u32 *)primes, count);
    memset(memory, 0, 1024);

    find_primes_sieve_seq(2, 500, primes);
    print_sieve(primes, 2, 500);
    memset(memory, 0, 1024);

    find_primes_sieve_seq_v2(2, 500, primes);
    print_sieve_v2(primes, 2, 500);
    memset(memory, 0, 1024);

    find_primes_sieve_seq_v3(2, 500, primes);
    print_sieve_v3(primes, 2, 500);
    memset(memory, 0, 1024);

    find_primes_sieve_seq_v4(2, 500, primes);
    print_sieve_v4(primes, 2, 500);
    memset(memory, 0, 1024);

    //find_primes_sieve_domain(2, 1000000000, primes);
    //find_primes_sieve_domain_v2(2000000000, 4000000000, primes);
    //find_primes_sieve_functional(2, 1000000000, primes);
    //find_primes_sieve_domain(2, 1000000000, primes);
    //find_primes_sieve_seq_v4(2, 4000000000, primes);
    return 0;
}
