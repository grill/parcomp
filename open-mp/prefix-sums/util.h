#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <assert.h>

#define ATYPE long int
#define ASIZE (1 << 20)
#define AVAL 1

#define STYPE long int
#define OTYPE long int

void hillis_steele_scan(ATYPE** x, STYPE n, OTYPE* ops, int threads);
void single_hillis_steele_scan(ATYPE** x, STYPE n, OTYPE* ops);
void hillis_steele_scan_time(ATYPE** x, STYPE n, int threads);
void single_hillis_steele_scan_time(ATYPE** x, STYPE n);

void single_recursive_scan(ATYPE** x, STYPE n, OTYPE* ops);
void recursive_scan(ATYPE** x, STYPE n, OTYPE* ops, int threads);
void single_recursive_scan_time(ATYPE** x, STYPE n);
void recursive_scan_time(ATYPE** x, STYPE n, int threads);
ATYPE recursive_scan_sum(ATYPE** x, STYPE n, OTYPE* ops, int threads);

void single_iterative_scan(ATYPE** x, STYPE n, OTYPE* ops);
void iterative_scan(ATYPE** x, STYPE n, OTYPE* ops, int threads);
void single_iterative_scan_time(ATYPE** x, STYPE n);
void iterative_scan_time(ATYPE** x, STYPE n, int threads);
ATYPE iterative_scan_sum(ATYPE** x, STYPE n, OTYPE* ops, int threads);