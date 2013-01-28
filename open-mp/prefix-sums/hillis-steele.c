#include "util.h"

// fast, not work-optimal (but still usefull)

//O ( n log n )

void hillis_steele_scan_time (ATYPE** x, STYPE n, int threads) {
  ATYPE *y = calloc(n, sizeof(ATYPE));
  ATYPE *t;
  STYPE i;
  STYPE k;

  for (k=1; k<n; k<<=1) {
    #pragma omp parallel num_threads(threads)
    {
       #pragma omp for nowait
       for (i=0; i<k; i++) {
         y[i] = (*x)[i];
       }
       #pragma omp for nowait
       for (i=k; i<n; i++) {
         y[i] = (*x)[i-k]+(*x)[i];
       }
    }
    t = (*x); (*x) = y; y = t; // swap
  }

  free(y);
}

void hillis_steele_scan(ATYPE** x, STYPE n, OTYPE* ops, int threads) {
  ATYPE *y = calloc(n, sizeof(ATYPE));
  ATYPE *t;
  STYPE i;
  STYPE k;

  for (k=1; k<n; k<<=1) {
    #pragma omp parallel num_threads(threads)
    {
       #pragma omp for nowait
       for (i=0; i<k; i++) {
         y[i] = (*x)[i];
         ops[omp_get_thread_num()]++;
       }
       #pragma omp for nowait
       for (i=k; i<n; i++) {
         y[i] = (*x)[i-k]+(*x)[i];
         ops[omp_get_thread_num()]++;
       }
    }
    t = (*x); (*x) = y; y = t; // swap
  }

  free(y);
}

void single_hillis_steele_scan_time (ATYPE** x, STYPE n) {
  ATYPE *y = calloc(n, sizeof(ATYPE));
  ATYPE *t;
  STYPE i;
  STYPE k;

  for (k=1; k<n; k<<=1) {
    for (i=0; i<k; i++) {
      y[i] = (*x)[i];
    }
    for (i=k; i<n; i++) {
      y[i] = (*x)[i-k]+(*x)[i];
    }
    t = (*x); (*x) = y; y = t; // swap
  }

  free(y);
}

void single_hillis_steele_scan(ATYPE** x, STYPE n, OTYPE* ops) {
  ATYPE *y = calloc(n, sizeof(ATYPE));
  ATYPE *t;
  STYPE i;
  STYPE k;

  for (k=1; k<n; k<<=1) {
    for (i=0; i<k; i++) {
      y[i] = (*x)[i];
      ops[omp_get_thread_num()]++;
    }
    for (i=k; i<n; i++) {
      y[i] = (*x)[i-k]+(*x)[i];
      ops[omp_get_thread_num()]++;
    }
    t = (*x); (*x) = y; y = t; // swap
  }

  free(y);
}