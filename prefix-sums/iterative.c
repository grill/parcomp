
#include <stdio.h>
#include <stdlib.h>
#include "util.h"

// OpenMP header
#include <omp.h>


int main(int argc, char *argv[]) {
  int threads, myid;
  int i;

  threads = 1;
  for (i=1; i<argc&&argv[i][0]=='-'; i++) {
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
