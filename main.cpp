#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <omp.h>

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
    u32 stub;
    u32 base_prime_count = find_primes_seq(2, sqrt_end, base_primes, &stub);

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


u32 find_primes_parallel(u32 begin, u32 end, u32 *primes) {
    u32 prime_count = 0;

    u32 sqrt_end = sqrt(end);
    u32 stub;
    prime_count = find_primes_seq(2, sqrt_end, primes, &stub);

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


    return prime_count;
}

u32 find_primes_parallel_v3(u32 begin, u32 end, u32 *primes) {
    u32 prime_count = 0;

    u32 sqrt_end = sqrt(end);
    u32 stub;
    prime_count = find_primes_seq(2, sqrt_end, primes, &stub);

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
    u32 stub;
    u32 base_prime_count = find_primes_seq(2, sqrt_end, base_primes, &stub);
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

#pragma omp parallel shared(primes)
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

#pragma omp parallel shared(primes, base_primes)
    {
        u32 work_per_thread = base_prime_count / omp_get_num_threads() + 1;
        u32 thread_begin = work_per_thread*omp_get_thread_num();
        u32 thread_end = thread_begin + work_per_thread;
        if (thread_end > base_prime_count)
            thread_end = base_prime_count;
#if 0
        printf("thread %d starts at %d(%d), ends at %d(%d)\n", omp_get_thread_num(), thread_begin, base_primes[thread_begin], thread_end-1, base_primes[thread_end-1]);
#endif

        for (u32 i = thread_begin; i < thread_end; i++) {
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

#define get_bit(a, i) (a[i>>3] & (1 << (i&7)))
#define set_bit(a, i) (a[i>>3] |= (1 << (i&7)))

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
            u32 start = begin-1 + i - ((begin-1) % i);
            start += i * ((start & 1) == 0);
            for (u32 j = start; j <= end; j += i+i) {
                set_bit(primes, j>>1);
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
    omp_set_num_threads(4);

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

#if 0
        for (u32 i = 3; i <= sqrt_end; i+=2) {
            if (!get_bit(primes, i>>1)) {
                u32 start = thread_begin-1 + i - ((thread_begin-1) % i);
                start += i * ((start & 1) == 0);
                for (u32 j = start; j < thread_end; j += i+i) {
                    set_bit(primes, j>>1);
                }
            }
        }
#else
        find_primes_sieve_seq_v4_no_init(thread_begin, thread_end-1, primes);
#endif
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
    test_sieve_domain_v2(memory);
#if 0
    auto start_time = omp_get_wtime();
    find_primes_sieve_domain_v2(2, 4000000000, primes);
    auto end_time = omp_get_wtime();
    u32 prime_count = 0;
    for (u32 i = 2; i <= 4000000000; i++) {
        if (i==2 || ((i&1) && !get_bit(primes, i/2))) {
            prime_count++;
        }
    }
    printf("found %d primes in %f seconds\n", prime_count, end_time-start_time);
#else
    //find_primes_sieve_domain_v2(2, 4000000000, primes);
    //find_primes_sieve_domain(2, 1000000000, primes);
    //find_primes_sieve_seq_v4(2, 4000000000, primes);
#endif

#if 0
#pragma omp parallel sections
    {
#pragma omp section
        find_primes_sieve_seq_v4(31616, 250023719, primes);
#pragma omp section
        find_primes_sieve_seq_v4(250023720, 500015823, primes);
#pragma omp section
        find_primes_sieve_seq_v4(500015824, 750007927, primes);
#pragma omp section
        find_primes_sieve_seq_v4(750007928, 1000000000, primes);
    }
#endif

    return 0;
}
