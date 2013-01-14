#include <stdio.h>
#include <stdlib.h>
#include "util.h"

// OpenMP header
#include <omp.h>

#define ATYPE int

int main(int argc, char *argv[]) {
  int threads, myid;
  int i;

  threads = 1;
  for (i = 1; i<argc && argv[i][0] == '-'; i++) {
    if (argv[i][1]=='t') i++,sscanf(argv[i],"%d",&threads);
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

  return 0;
}

void Scan(ATYPE[] x, int n) {
  if(n==1) return;

  ATYPE* y = malloc(sizeof(ATYPE)*(n/2));

  #pragma omp parallel for
  for(i=0; i<n/2; i++) y[i] = x[2*i]+x[2*i+1];

  hScan(y, n/2);

  x[1] = y[0];
  #pragma omp parallel for
  for(int i = 1; i < n/2; i++) {
    x[2*i] = y[i-1]+x[2*i];
    x[2*i+1] = y[i];
  }

  if (odd(n)) x[n-1] = y[n/2-1]+x[n-1];

  free(y);
}