#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "util.h"

// OpenMP header
#include <omp.h>

#define ATYPE int
#define ASIZE 100000000
#define AVAL 1

int* gen_int_array(int v, int n);


long big_sum(int* x, int n);
int reduce_scan(ATYPE* x, int n);
int single_sum(ATYPE* x, int n);

void single_recursive_scan(ATYPE** x, int n, int* ops);
void recursive_scan(ATYPE** x, int n, int* ops);
void iterative_scan(ATYPE** x, int n, int* ops);
void single_iterative_scan(ATYPE** x, int n, int* ops);
void hillis_steele_scan(ATYPE** x, int n, int* ops);
void single_hillis_steele_scan(ATYPE** x, int n, int* ops);

void testReduce(char* name, int (*)(ATYPE*, int));
void testPrefixsum(char* name, void (*)(ATYPE**, int, int*));

//void powers_iterative_scan(ATYPE** x, int n);

int main(int argc, char *argv[]) {
 // int threads, myid;
 // int i;
 /*
  int* arr_reduc;
  int* arr_recur;
  int* arr_single;
  int* arr_single_sum;
  int sum_reduc, sum_recur, sum_single, sum_single_sum;

  double time_reduc_1, time_reduc_2, time_recur_1, time_recur_2;
  double time_single_1, time_single_2, time_single_sum_1, time_single_sum_2;
*/
 /*
  threads = 1;
  for (i = 1; i<argc && argv[i][0] == '-'; i++) {
    if (argv[i][1]=='t') i++, sscanf(argv[i],"%d",&threads);
  }

  printf("Maximum number of threads possible is %d\n",omp_get_max_threads());

  if (threads<omp_get_max_threads()) {
    if (threads<1) threads = 1;
    omp_set_num_threads(threads);
  } else {
    threads = omp_get_max_threads();
  }

  #pragma omp parallel num_threads(threads)
  //#pragma omp single
  {
    myid = omp_get_thread_num();
    printf("Thread %d of %d active\n",myid,threads);
  }
  */
/*
  arr_reduc = gen_int_array(1, ASIZE);

  time_reduc_1 = omp_get_wtime();
  sum_reduc = reduce_scan(arr_reduc, ASIZE);
  time_reduc_2 = omp_get_wtime();

  arr_recur = gen_int_array(1, ASIZE);

  time_recur_1 = omp_get_wtime();
  recursive_scan(&arr_recur, ASIZE);
  time_recur_2 = omp_get_wtime();
  sum_recur = arr_recur[ASIZE-1];

  arr_single = gen_int_array(1, ASIZE);

  time_single_1 = omp_get_wtime();
  recursive_scan(&arr_single, ASIZE);
  time_single_2 = omp_get_wtime();
  sum_single = arr_single[ASIZE-1];

  arr_single_sum = gen_int_array(1, ASIZE);

  time_single_sum_1 = omp_get_wtime();
  sum_single_sum = single_sum(arr_single_sum, ASIZE);
  time_single_sum_2 = omp_get_wtime();

  fprintf(stdout, "Parallel Reduction -- startTime: %f, endTime: %f, diff: %f, sum: %i\n",
             time_reduc_1, time_reduc_2, (time_reduc_2-time_reduc_1), sum_reduc);

  fprintf(stdout, "Single Reduction -- startTime: %f, endTime: %f, diff: %f, sum: %i\n",
             time_single_sum_1, time_single_sum_2, (time_single_sum_2-time_single_sum_1), sum_single_sum);

  fprintf(stdout, "Parallel Recursion -- startTime: %f, endTime: %f, diff: %f, sum: %i\n",
             time_recur_1, time_recur_2, (time_recur_2-time_recur_1), sum_recur);

  fprintf(stdout, "Single Recursion -- startTime: %f, endTime: %f, diff: %f, sum: %i\n",
             time_single_1, time_single_2, (time_single_2-time_single_1), sum_single);
  */

  //anzahl befehle
  //speed with threadnum

  //small data
  //big data

  

  //unterschiedliche Daten ??

  //-t threads a von null rauf zaehlen -p wie oft durcheglaufen wird
  //calloc statt malloc
  //TODO: Clean Up, Split Count ops & performance, hills-steele performance after malloc

  double size = ASIZE;

  fprintf(stdout, "--- Stats ---\n\n");
  fprintf(stdout, "n: %f, n*log(n): %f, n*log(n)^2: %f, n*n: %f\n\n",
	  size, size*log2(size), size*log2(size)*log2(size), size*size);

  testReduce("Parallel Reduction", &reduce_scan);
  testReduce("Single Reduction", &single_sum);
  fprintf(stdout, "\n");

  testPrefixsum("Parallel Recursion", &recursive_scan);
  testPrefixsum("Single Recursion", &single_recursive_scan);
  fprintf(stdout, "\n");

  testPrefixsum("Parallel Iterative", &iterative_scan);
  //testPrefixsum("Parallel Iterative with powers of 2", &powers_iterative_scan);
  testPrefixsum("Single Iterative", &single_iterative_scan);
  fprintf(stdout, "\n");

  testPrefixsum("Parallel Hillis-Steele", &hillis_steele_scan);
  testPrefixsum("Single Hillis-Steele", &single_hillis_steele_scan);

  return 0;
}

void testReduce(char* name, int (*sum)(ATYPE*, int)) {
  ATYPE* arr;
  double time_begin, time_end;
  int erg;

  arr = gen_int_array(AVAL, ASIZE);

  time_begin = omp_get_wtime();
  erg = sum(arr, ASIZE);
  time_end = omp_get_wtime();

  fprintf(stdout, "%s -- startTime: %f, endTime: %f, diff: %f, sum: %i\n",
             name, time_begin, time_end, (time_end - time_begin), erg);
}

void testPrefixsum(char* name, void (*prefixsum)(ATYPE**, int, int*)) {
  ATYPE* arr;
  double time_begin, time_end;
  int* ops = malloc(sizeof(int) * omp_get_max_threads());
  memset((void *) ops, 0, sizeof(int)*omp_get_max_threads());

  arr = gen_int_array(AVAL, ASIZE);

  time_begin = omp_get_wtime();
  prefixsum(&arr, ASIZE, ops);
  time_end = omp_get_wtime();

  fprintf(stdout, "%s -- ops: %li startTime: %f, endTime: %f, diff: %f, sum: %i\n",
          name, big_sum(ops, omp_get_max_threads()),
          time_begin, time_end, (time_end - time_begin), arr[ASIZE-1]);

  free(ops);
  free(arr);
}

int* gen_int_array(int v, int n) {
  int* ret = malloc(sizeof(int)*n);

  for(int i = 0; i < n; i++)
    if(i % 2 == 0)
       ret[i] = v;

  return ret;
}


int single_sum(ATYPE* x, int n) {
  int sum = 0;

  for(int i = 0; i < n;i++) {
     sum += x[i];
  }

  return sum;
}

long big_sum(ATYPE* x, int n) {
  long sum = 0;

  #pragma omp parallel for reduction(+:sum)
  for(int i = 0; i < n;i++) {
     sum += x[i];
  }

  return sum;
}

int reduce_scan(ATYPE* x, int n) {
  int sum = 0;

  #pragma omp parallel for reduction(+:sum)
  for(int i = 0; i < n;i++) {
     sum += x[i];
  }

  return sum;
}

void single_recursive_scan(ATYPE** x, int n, int* ops) {
  int i;
  if(n==1) return;

  ATYPE* y = malloc(sizeof(ATYPE)*(n/2));

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

void recursive_scan(ATYPE** x, int n, int* ops) {
  int i;

  if(n==1) return;

  ATYPE* y = malloc(sizeof(ATYPE)*(n/2));

  #pragma omp parallel for
  for(i=0; i<n/2; i++) {
    y[i] = (*x)[2*i]+(*x)[2*i+1];
    ops[omp_get_thread_num()]++;
  }

  recursive_scan(&y, n/2, ops);

  (*x)[1] = y[0];
  #pragma omp parallel for
  for(i = 1; i < n/2; i++) {
   (*x)[2*i] = y[i-1]+(*x)[2*i];
   (*x)[2*i+1] = y[i];
    ops[omp_get_thread_num()]+=2;
  }

  if (n % 2 != 0) {
    (*x)[n-1] = y[n/2-1]+(*x)[n-1];
    ops[omp_get_thread_num()]++;
  }

  free(y);
}

/*
void powers_iterative_scan(ATYPE** x, int n) {
  int k=0;
  int kk;
  int i;
  int max = log2(n);
  int k2;

  //up phase
  #pragma omp parallel for
  for (k=0; k<max; k++) {
    k2 = pow(2, k);
    kk = k2<<1; // double
    for (i=kk-1; i<n; i+=kk) {
       (*x)[i] = (*x)[i-k2]+(*x)[i];
    }
  }

  //down phase
  #pragma omp parallel for
  for (k = k-1; k>1; k--) {
    k2 = pow(2, k);
    kk = k2>>1; // double
    for (i=k-1; i<n-kk; i+=k2) {
       (*x)[i+kk] = (*x)[i]+(*x)[i+kk];
    }
  }

}*/

void iterative_scan(ATYPE** x, int n, int* ops) {
  int k;
  int kk;
  int i;

  //up phase
  for (k=1; k<n; k=kk) {
    kk = k<<1; // double
    #pragma omp parallel for
    for (i=kk-1; i<n; i+=kk) {
       (*x)[i] = (*x)[i-k]+(*x)[i];
       ops[omp_get_thread_num()]++;
    }
  }

//down phase
  for (k=k>>1; k>1; k=kk) {
    kk = k>>1; // halve
    #pragma omp parallel for
    for (i=k-1; i<n-kk; i+=k) {
       (*x)[i+kk] = (*x)[i]+(*x)[i+kk];
       ops[omp_get_thread_num()]++;
    }
  }
}

void single_iterative_scan(ATYPE** x, int n, int* ops) {
  int k;
  int kk;
  int i;

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

//besser aufteilen
void hillis_steele_scan(ATYPE** x, int n, int* ops) {
  int *y = malloc(n*sizeof(int));
  int *t;
  int i;
  int k;

  for (k=1; k<n; k<<=1) {
    #pragma omp parallel for
    for (i=0; i<k; i++) {
      y[i] = (*x)[i];
      ops[omp_get_thread_num()]++;
    }
    #pragma omp parallel for
    for (i=k; i<n; i++) {
      y[i] = (*x)[i-k]+(*x)[i];
      ops[omp_get_thread_num()]++;
    }
    t = (*x); (*x) = y; y = t; // swap
  }

  free(y);
}

void single_hillis_steele_scan(ATYPE** x, int n, int* ops) {
  int *y = malloc(n*sizeof(int));
  int *t;
  int i;
  int k;

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
