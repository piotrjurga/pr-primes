void test_find_seq(void *memory) {
    u32 *primes = (u32 *)memory;

    u32 begin;
    u32 count = find_primes_seq(2, 50, primes, &begin);
    puts("test find_primes_seq from 2 to 50");
    printf("found %d + %d\n", begin, count-begin);
    for (u32 i = 0; i < count; i++) {
        printf("%d;", primes[i]);
    }
    puts("");

    count = find_primes_seq(25, 50, primes, &begin);
    puts("test find_primes_seq from 25 to 50");
    printf("found %d + %d\n", begin, count-begin);
    for (u32 i = 0; i < count; i++) {
        printf("%d;", primes[i]);
    }
    puts("");

    auto start_time = omp_get_wtime();
    find_primes_seq(2, 10000000, primes, &begin);
    auto end_time = omp_get_wtime();
    printf("time: %f\n", end_time-start_time);
}

void test_find_seq_v2(void *memory) {
    u8 *primes = (u8 *)memory;

    find_primes_seq_v2(2, 50, primes);
    puts("test find_primes_seq_v2 from 2 to 50");
    for (u32 i = 2; i <= 50; i++) {
        if (!primes[i]) {
            printf("%d;", i);
        }
    }
    puts("");

    find_primes_seq_v2(25, 50, primes);
    puts("test find_primes_seq_v2 from 25 to 50");
    for (u32 i = 25; i <= 50; i++) {
        if (!primes[i]) {
            printf("%d;", i);
        }
    }
    puts("");

    auto start_time = omp_get_wtime();
    find_primes_seq_v2(2, 10000000, primes);
    auto end_time = omp_get_wtime();
    printf("time: %f\n", end_time-start_time);
}

void test_find_parallel(void *memory) {
    u32 *primes = (u32 *)memory;

    puts("test find_primes_parallel from 2 to 10000000");
    auto start_time = omp_get_wtime();
    u32 count = find_primes_parallel(2, 10000000, primes);
    auto end_time = omp_get_wtime();
    printf("%d\n", count);
    printf("time: %f\n", end_time-start_time);
/*
    for (u32 i = 0; i < count; i++) {
        printf("%d;", primes[i]);
    }
    puts("");
*/
}

void test_find_parallel_v2(void *memory) {
    u8 *primes = (u8 *)memory;

    puts("test find_primes_parallel_v2 from 2 to 50");
    find_primes_parallel_v2(2, 50, primes);
    for (u32 i = 2; i <= 50; i++) {
        if (!primes[i]) {
            printf("%d;", i);
        }
    }
    puts("");

    puts("test find_primes_parallel_v2 from 25 to 50");
    find_primes_parallel_v2(25, 50, primes);
    for (u32 i = 25; i <= 50; i++) {
        if (!primes[i]) {
            printf("%d;", i);
        }
    }
    puts("");

    puts("test find_primes_parallel_v2 from 2 to 10000000");
    auto start_time = omp_get_wtime();
    find_primes_parallel_v2(2, 10000000, primes);
    auto end_time = omp_get_wtime();
    u32 prime_count = 0;
    for (u32 i = 2; i <= 10000000; i++) {
        if (!primes[i]) {
            prime_count++;
        }
    }
    printf("found %d primes in %f seconds\n", prime_count, end_time-start_time);
    puts("");
}

void test_sieve_seq(void *memory) {
    u8 *primes = (u8 *)memory;

    puts("test find_primes_sieve_seq from 2 to 50");
    find_primes_sieve_seq(2, 50, primes);
    for (u32 i = 2; i <= 50; i++) {
        if (!primes[i]) {
            printf("%d;", i);
        }
    }
    puts("");

    memset(primes, 0, 50);
    puts("test find_primes_sieve_seq from 25 to 50");
    find_primes_sieve_seq(25, 50, primes);
    for (u32 i = 25; i <= 50; i++) {
        if (!primes[i]) {
            printf("%d;", i);
        }
    }
    puts("");
}

void test_sieve_domain(void *memory) {
    u8 *primes = (u8 *)memory;

    puts("test find_primes_sieve_domain from 2 to 50");
    find_primes_sieve_domain(2, 50, primes);
    for (u32 i = 2; i <= 50; i++) {
        if (!primes[i]) {
            printf("%d;", i);
        }
    }
    puts("");

    memset(primes, 0, 50);
    puts("test find_primes_sieve_domain from 25 to 50");
    find_primes_sieve_domain(25, 50, primes);
    for (u32 i = 25; i <= 50; i++) {
        if (!primes[i]) {
            printf("%d;", i);
        }
    }
    puts("");

    puts("testing find_primes_sieve_domain from 2 to 1000000000");
    memset(primes, 0, 50);
    auto start_time = omp_get_wtime();
    find_primes_sieve_domain(2, 1000000000, primes);
    auto end_time = omp_get_wtime();
    u32 prime_count = 0;
    for (u32 i = 2; i <= 1000000000; i++) {
        if (!primes[i]) prime_count++;
    }
    printf("found %d primes in %f seconds\n", prime_count, end_time-start_time);
}

void test_sieve_functional(void *memory) {
    u8 *primes = (u8 *)memory;

    puts("test find_primes_sieve_functional from 2 to 50");
    find_primes_sieve_functional(2, 50, primes);
    for (u32 i = 2; i <= 50; i++) {
        if (!primes[i]) {
            printf("%d;", i);
        }
    }
    puts("");

    memset(primes, 0, 50);
    puts("test find_primes_sieve_functional from 25 to 50");
    find_primes_sieve_functional(25, 50, primes);
    for (u32 i = 25; i <= 50; i++) {
        if (!primes[i]) {
            printf("%d;", i);
        }
    }
    puts("");

    puts("testing find_primes_sieve_functional from 2 to 1000000000");
    memset(primes, 0, 1024*1024*1024);
    auto start_time = omp_get_wtime();
    find_primes_sieve_functional(2, 1000000000, primes);
    auto end_time = omp_get_wtime();
    u32 prime_count = 0;
    for (u32 i = 2; i <= 1000000000; i++) {
        if (!primes[i]) prime_count++;
    }
    printf("found %d primes in %f seconds\n", prime_count, end_time-start_time);
}


void test_sieve_seq_v2(void *memory) {
    u8 *primes = (u8 *)memory;

    puts("test find_primes_sieve_seq_v2 from 2 to 50");
    find_primes_sieve_seq_v2(2, 50, primes);
    for (u32 i = 2; i <= 50; i++) {
        if ((i&1) && !primes[i/2]) {
            printf("%d;", i);
        }
    }
    puts("");

    memset(primes, 0, 50);
    puts("test find_primes_sieve_seq_v2 from 25 to 50");
    find_primes_sieve_seq_v2(25, 50, primes);
    for (u32 i = 25; i <= 50; i++) {
        if ((i&1) && !primes[i/2]) {
            printf("%d;", i);
        }
    }
    puts("");
}

void test_sieve_seq_v3(void *memory) {
    u8 *primes = (u8 *)memory;

    puts("test find_primes_sieve_seq_v3 from 2 to 50");
    find_primes_sieve_seq_v3(2, 50, primes);
    for (u32 i = 2; i <= 50; i++) {
        if (!get_bit(primes, i)) {
            printf("%d;", i);
        }
    }
    puts("");

    memset(primes, 0, 50);
    puts("test find_primes_sieve_seq_v3 from 25 to 50");
    find_primes_sieve_seq_v3(25, 50, primes);
    for (u32 i = 25; i <= 50; i++) {
        if (!get_bit(primes, i)) {
            printf("%d;", i);
        }
    }
    puts("");
}

void test_sieve_seq_v4(void *memory) {
    u8 *primes = (u8 *)memory;

    puts("test find_primes_sieve_seq_v4 from 2 to 50");
    find_primes_sieve_seq_v4(2, 50, primes);
    for (u32 i = 2; i <= 50; i++) {
        if ((i&1) && !get_bit(primes, i/2)) {
            printf("%d;", i);
        }
    }
    puts("");

    memset(primes, 0, 50);
    puts("test find_primes_sieve_seq_v4 from 25 to 50");
    find_primes_sieve_seq_v4(25, 50, primes);
    for (u32 i = 25; i <= 50; i++) {
        if ((i&1) && !get_bit(primes, i/2)) {
            printf("%d;", i);
        }
    }
    puts("");

    puts("test find_primes_sieve_seq_v4 from 2 to 1000000000");
    memset(primes, 0, 1000000000);
    auto start_time = omp_get_wtime();
    find_primes_sieve_seq_v4(2, 1000000000, primes);
    auto end_time = omp_get_wtime();
    u32 prime_count = 1;
    for (u32 i = 3; i < 1000000000; i++) {
        if ((i&1) && !get_bit(primes, i/2)) {
            prime_count++;
        }
    }
    printf("found %d primes in %f seconds\n", prime_count, end_time-start_time);

    puts("test find_primes_sieve_seq_v4 from 500000000 to 1000000000");
    memset(primes, 0, 1000000000);
    start_time = omp_get_wtime();
    find_primes_sieve_seq_v4(500000000, 1000000000, primes);
    end_time = omp_get_wtime();
    prime_count = 0;
    for (u32 i = 500000000; i <= 1000000000; i++) {
        if ((i&1) && !get_bit(primes, i/2)) {
            prime_count++;
        }
    }
    printf("found %d primes in %f seconds\n", prime_count, end_time-start_time);
}

void test_sieve_seq_v(void *memory) {
    u8 *primes = (u8 *)memory;

    memset(primes, 0, 100000000);
    auto start_time = omp_get_wtime();
    find_primes_sieve_seq(2, 100000000, primes);
    auto end_time = omp_get_wtime();
    printf("v1 time = %f\n", end_time-start_time);

    memset(primes, 0, 100000000);
    start_time = omp_get_wtime();
    find_primes_sieve_seq_v2(2, 100000000, primes);
    end_time = omp_get_wtime();
    printf("v2 time = %f\n", end_time-start_time);

    memset(primes, 0, 100000000);
    start_time = omp_get_wtime();
    find_primes_sieve_seq_v3(2, 100000000, primes);
    end_time = omp_get_wtime();
    printf("v3 time = %f\n", end_time-start_time);

    memset(primes, 0, 100000000);
    start_time = omp_get_wtime();
    find_primes_sieve_seq_v4(2, 100000000, primes);
    end_time = omp_get_wtime();
    printf("v4 time = %f\n", end_time-start_time);
}

void test_sieve_domain_v2(void *memory) {
    u8 *primes = (u8 *)memory;

    memset(primes, 0, 50);
    puts("test find_primes_sieve_domain_v2 from 2 to 50");
    find_primes_sieve_domain_v2(2, 50, primes);
    for (u32 i = 2; i <= 50; i++) {
        if (i==2 || ((i&1) && !get_bit(primes, i/2))) {
            printf("%d;", i);
        }
    }
    puts("");

    memset(primes, 0, 50);
    puts("test find_primes_sieve_domain_v2 from 25 to 50");
    find_primes_sieve_domain_v2(25, 50, primes);
    for (u32 i = 25; i <= 50; i++) {
        if ((i&1) && !get_bit(primes, i/2)) {
            printf("%d;", i);
        }
    }
    puts("");

    u32 max = 100000000;

    puts("testing find_primes_sieve_domain_v2 from 2 to 1000000000");
    memset(primes, 0, 1000000000);
    auto start_time = omp_get_wtime();
    find_primes_sieve_domain_v2(2, max, primes);
    auto end_time = omp_get_wtime();
    u32 prime_count = 0;
    for (u32 i = 2; i <= max; i++) {
        if (i==2 || ((i&1) && !get_bit(primes, i/2))) {
            prime_count++;
        }
    }
    printf("found %d primes in %f seconds\n", prime_count, end_time-start_time);

    puts("testing find_primes_sieve_domain_v2 from 2 to 1000000000");
    memset(primes, 0, 1000000000);
    start_time = omp_get_wtime();
    find_primes_sieve_seq_v4(2, max, primes);
    end_time = omp_get_wtime();
    prime_count = 0;
    for (u32 i = 2; i <= max; i++) {
        if (i==2 || ((i&1) && !get_bit(primes, i/2))) {
            prime_count++;
        }
    }
    printf("found %d primes in %f seconds\n", prime_count, end_time-start_time);


#if 0
    memset(primes, 0, 1000000000);
    start_time = omp_get_wtime();
    find_primes_sieve_seq_v4(2, 1000000000, primes);
    end_time = omp_get_wtime();
    prime_count = 0;
    for (u32 i = 2; i <= 1000000000; i++) {
        if (i==2 || ((i&1) && !get_bit(primes, i/2))) {
            prime_count++;
        }
    }
    printf("found %d primes in %f seconds\n", prime_count, end_time-start_time);
#endif
#if 0
    puts("testing hardcoded parallism from 2 to 1000000000");
    memset(primes, 0, 1000000000);
    start_time = omp_get_wtime();
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
    end_time = omp_get_wtime();
    prime_count = 0;
    for (u32 i = 2; i <= 1000000000; i++) {
        if (i==2 || ((i&1) && !get_bit(primes, i/2))) {
            prime_count++;
        }
    }
    printf("found %d primes in %f seconds\n", prime_count, end_time-start_time);
#endif
}
