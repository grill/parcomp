#include "util.h"

// fast, work-optimal

// T(n) = log n
// T fast (n) = log n
// O (n)
// T seq (n) = O( n )
// T par (p, n) = O ( n / p ) + T fast (n)
// Linear speed-up up to T seq (n) / T fast (n) processors, n / log n

// Drawbacks:
// Space: extra n / 2 sized array for each recursive call, n in total
// About 2n opterations (seq scan: n-1)
// 2 log2 n parallel steps

void single_recursive_scan_time (ATYPE** x, STYPE n) {
  STYPE i;
  if(n==1) return;

  ATYPE* y = calloc(n/2, sizeof(ATYPE));

  for(i=0; i<n/2; i++) { 
    y[i] = (*x)[2*i]+(*x)[2*i+1];
  }

  single_recursive_scan_time(&y, n/2);

  (*x)[1] = y[0];
  for(i = 1; i < n/2; i++) {
    (*x)[2*i] = y[i-1]+(*x)[2*i];
    (*x)[2*i+1] = y[i];
  }

  if (n % 2 != 0) {
    (*x)[n-1] = y[n/2-1]+(*x)[n-1];
  }

  free(y);
}

void single_recursive_scan(ATYPE** x, STYPE n, OTYPE* ops) {
  STYPE i;
  if(n==1) return;

  ATYPE* y = calloc(n/2, sizeof(ATYPE));

  for(i=0; i<n/2; i++) { 
    y[i] = (*x)[2*i]+(*x)[2*i+1];
    ops[omp_get_thread_num()]++;
  }

  single_recursive_scan(&y, n/2, ops);

  (*x)[1] = y[0];
  for(i = 1; i < n/2; i++) {
    (*x)[2*i] = y[i-1]+(*x)[2*i];
    (*x)[2*i+1] = y[i];
    ops[omp_get_thread_num()] += 2;
  }

  if (n % 2 != 0) {
    (*x)[n-1] = y[n/2-1]+(*x)[n-1];
    ops[omp_get_thread_num()]++;
  }

  free(y);
}

void recursive_scan_time (ATYPE** x, STYPE n, int threads) {
  STYPE i;

  if(n==1) return;

  ATYPE* y = calloc(n/2, sizeof(ATYPE));

  #pragma omp parallel for num_threads(threads)
  for(i=0; i<n/2; i++) {
    y[i] = (*x)[2*i]+(*x)[2*i+1];
  }

  recursive_scan_time(&y, n/2, threads);

  (*x)[1] = y[0];
  #pragma omp parallel for num_threads(threads)
  for(i = 1; i < n/2; i++) {
   (*x)[2*i] = y[i-1]+(*x)[2*i];
   (*x)[2*i+1] = y[i];
  }

  if (n % 2 != 0) {
    (*x)[n-1] = y[n/2-1]+(*x)[n-1];
  }

  free(y);
}

void recursive_scan(ATYPE** x, STYPE n, OTYPE* ops, int threads) {
  STYPE i;

  if(n==1) return;

  ATYPE* y = calloc(n/2, sizeof(ATYPE));

  #pragma omp parallel for num_threads(threads)
  for(i=0; i<n/2; i++) {
    y[i] = (*x)[2*i]+(*x)[2*i+1];
    ops[omp_get_thread_num()]+=2;
  }

  recursive_scan(&y, n/2, ops, threads);

  (*x)[1] = y[0];
  #pragma omp parallel for num_threads(threads)
  for(i = 1; i < n/2; i++) {
   (*x)[2*i] = y[i-1]+(*x)[2*i];
   (*x)[2*i+1] = y[i];
    ops[omp_get_thread_num()]+=3;
  }

  if (n % 2 != 0) {
    (*x)[n-1] = y[n/2-1]+(*x)[n-1];
    ops[omp_get_thread_num()]++;
  }

  free(y);
}

ATYPE recursive_scan_sum(ATYPE** x, STYPE n, OTYPE* ops, int threads) {
  STYPE i;

  if(n==1) return (*x)[0];

  ATYPE* y = calloc(n/2, sizeof(ATYPE));

  #pragma omp parallel for num_threads(threads)
  for(i=0; i<n/2; i++) {
    y[i] = (*x)[2*i]+(*x)[2*i+1];
    ops[omp_get_thread_num()]+=2;
  }

  recursive_scan_sum(&y, n/2, ops, threads);

  if((n % 2)  ^ 0) {
    (*x)[n-1] += y[n/2-1];
    ops[0]+=2;
  } else {
    (*x)[n-1] = y[n/2-1];
    ops[0]++;
  }

  free(y);
  return (*x)[n-1];
}