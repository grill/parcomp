#include "util.h"

// fast, work-optimal

// Total work ca. 2n = O ( T seq (n) )
// Speedup (p) at most p / 2
// For p=n: work optimal (not cost optimal) - p processors in 2 log p rounds
// = O (p log p)

void iterative_scan_time (ATYPE** x, STYPE n, int threads) {
  STYPE k;
  STYPE kk;
  STYPE i;

  //up phase
  for (k=1; k<n; k=kk) {
    kk = k<<1; // double
    #pragma omp parallel for num_threads(threads)
    for (i=kk-1; i<n; i+=kk) {
       (*x)[i] = (*x)[i-k]+(*x)[i];
    }
  }

//down phase
  for (k=k>>1; k>1; k=kk) {
    kk = k>>1; // halve
    #pragma omp parallel for num_threads(threads)
    for (i=k-1; i<n-kk; i+=k) {
       (*x)[i+kk] = (*x)[i]+(*x)[i+kk];
    }
  }
}

void iterative_scan(ATYPE** x, STYPE n, OTYPE* ops, int threads) {
  STYPE k;
  STYPE kk;
  STYPE i;

  //up phase
  for (k=1; k<n; k=kk) {
    kk = k<<1; // double
    #pragma omp parallel for num_threads(threads)
    for (i=kk-1; i<n; i+=kk) {
       (*x)[i] = (*x)[i-k]+(*x)[i];
       ops[omp_get_thread_num()]++;
    }
  }

//down phase
  for (k=k>>1; k>1; k=kk) {
    kk = k>>1; // halve
    #pragma omp parallel for num_threads(threads)
    for (i=k-1; i<n-kk; i+=k) {
       (*x)[i+kk] = (*x)[i]+(*x)[i+kk];
       ops[omp_get_thread_num()]++;
    }
  }
}

//Bonus: currently works only when n = 2^k
ATYPE iterative_scan_sum(ATYPE** x, STYPE n, OTYPE* ops, int threads) {
  STYPE k;
  STYPE kk;
  STYPE i;

  //up phase
  for (k=1; k<n; k=kk) {
    kk = k<<1; // double
    #pragma omp parallel for num_threads(threads)
    for (i=kk-1; i<n; i+=kk) {
       (*x)[i] = (*x)[i-k]+(*x)[i];
       ops[omp_get_thread_num()]++;
    }
  }

  if(k ^ n)
    for (k=k>>1,kk=k; kk<n && k>1;kk+=k) {
      (*x)[n-1] += (*x)[kk-1];
      ops[0]++;
      for(;k>1 && kk+k > n;k=k>>1);
    }
/*
  for (k=k>>1; k>1; k=kk) {
    kk = k>>1; // halve
    printf("k: %li, kk: %li, x[(k-1)]: %li, x[(k-1)+kk]: %li, idx1: %li, idx2: %li happens: %d\n",
	   k, kk, (*x)[(k-1)], (*x)[(k-1)+kk], (k-1)+kk, (k-1), k-1<n-kk);
    if(k-1<n-kk) {
      (*x)[(k-1)+kk] = (*x)[(k-1)]+(*x)[(k-1)+kk];
    }
  }*/

  return (*x)[n-1];
/*
  if(n - k == 0) {
    return (*x)[n-1];
  } else {
    k = k>>1;
    (*x) = (*x)+(k-1);
    printf("\nn: %li, k: %li, d: %li, d2: %li\n", n,k,(*x)[0], (*x)[1]);
    ret =  iterative_scan_sum(x, n-(k-1), ops, threads);
    (*x) = (*x)-(k-1);
    return ret;
  }*/
}

void single_iterative_scan_time (ATYPE** x, STYPE n) {
  STYPE k;
  STYPE kk;
  STYPE i;

  //up phase
  for (k=1; k<n; k=kk) {
    kk = k<<1; // double
    for (i=kk-1; i<n; i+=kk) {
       (*x)[i] = (*x)[i-k]+(*x)[i];
    }
  }

//down phase
  for (k=k>>1; k>1; k=kk) {
    kk = k>>1; // halve
    for (i=k-1; i<n-kk; i+=k) {
       (*x)[i+kk] = (*x)[i]+(*x)[i+kk];
    }
  }
}

void single_iterative_scan(ATYPE** x, STYPE n, OTYPE* ops) {
  STYPE k;
  STYPE kk;
  STYPE i;

  //up phase
  for (k=1; k<n; k=kk) {
    kk = k<<1; // double
    for (i=kk-1; i<n; i+=kk) {
       (*x)[i] = (*x)[i-k]+(*x)[i];
       ops[omp_get_thread_num()]++;
    }
  }

//down phase
  for (k=k>>1; k>1; k=kk) {
    kk = k>>1; // halve
    for (i=k-1; i<n-kk; i+=k) {
       (*x)[i+kk] = (*x)[i]+(*x)[i+kk];
       ops[omp_get_thread_num()]++;
    }
  }
}
