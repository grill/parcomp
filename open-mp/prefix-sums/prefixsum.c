#include "util.h"

ATYPE* gen_int_array(ATYPE v, STYPE n);
ATYPE big_sum(ATYPE* x, STYPE n);
ATYPE reduce_scan(ATYPE* x, STYPE n, int threads);
ATYPE single_sum(ATYPE* x, STYPE n);

void testReduce(char* name, STYPE n, ATYPE value, ATYPE erg, int correct, ATYPE (*)(ATYPE*, STYPE));
void testParaReduce(char* name, STYPE n, ATYPE value, ATYPE erg, int correct, int threads, ATYPE (*sum)(ATYPE*, STYPE, int));
void testPrefixsum(char* name, STYPE n, ATYPE value, int correct, void (*)(ATYPE**, STYPE, OTYPE*));
void testParaPrefixsum(char* name, STYPE n, ATYPE value, int correct ,int threads, void (*)(ATYPE**, STYPE, OTYPE*, int));
void testParaSum(char* name, STYPE n, ATYPE value, ATYPE sum, int correct, int threads, ATYPE (*prefixsum)(ATYPE**, STYPE, OTYPE*, int));
void testPrefixsum_time (char* name, STYPE n, ATYPE value, int correct, void (*prefixsum)(ATYPE**, STYPE));
void testParaPrefixsum_time (char* name, STYPE n, ATYPE value, int correct, int threads, void (*prefixsum)(ATYPE**, STYPE, int));

static void usage(void);
void printArray(ATYPE * arr, ATYPE size, ATYPE i);

int check_correctness(ATYPE* x, STYPE n);

  //anzahl befehle
  //speed with threadnum

  //small data
  //big data

  //TODO: Correctness by Testing missing; Bonus not fully fledged now (n^2 and some ifs); never crash
int main(int argc, char *argv[]) {
  int i;
  int threads = 1;
  int single = 0;
  int max_threads = 0;
  STYPE n = ASIZE;
  ATYPE value = AVAL;
  ATYPE* arr;

  int all_threads = 0;
  int correctness = 0;
  int performance = 0;
  int work = 0;

  int recursive = 0;
  int iterative = 0;
  int hills = 0;
  int reduction = 0;
  int only_sum = 0;

  ATYPE erg = 0;

  char c;
  char* tailptr;

   opterr = 0;
   while ((c = getopt(argc, argv, "n:v:t:smacpwriheo")) != -1 ) {
      switch (c) {
      case 'n': //array length
         if(sscanf(optarg,"%ld",&n) == EOF) {
            (void) fprintf(stderr, "%s: %s is not a valid number.\n", argv[0],
               optarg);
            usage();
         }
         break;
      case 'v': //value
         if(sscanf(optarg,"%ld",&value) == EOF) {
            (void) fprintf(stderr, "%s: %s is not a valid number.\n", argv[0],
               optarg);
            usage();
         }
         break;

      case 't': //thread count
         threads = strtol(optarg, &tailptr, 10);
         if(*tailptr != '\0') {
            (void) fprintf(stderr, "%s: %s is not a valid number.\n", argv[0],
               tailptr);
            usage();
         }
         break;
      case 's': /* single? */ single = 1; break;
      case 'm': /* max thread num */ max_threads = 1; break;
      case 'a': /* all threads */ all_threads = 1; break;

      case 'c': /* Corectness check */ correctness = 1; break;
      case 'p': /* Performance check */ performance = 1; break;
      case 'w': /* work check */ work = 1; break;

      case 'r': recursive = 1; break;
      case 'i': iterative = 1; break;
      case 'h': hills = 1; break;
      case 'e': reduction = 1; break;
      case 'o': only_sum = 1; break;

      case '?':
         (void) fprintf(stderr, "%s: The Parameter %c is invalid.\n", argv[0],
            optopt);
         usage();
         break;
      default:
         assert(0);
         break;
      }
   }

  if(all_threads || max_threads || threads > omp_get_max_threads()) {
     threads = omp_get_max_threads();
  } else if (threads < omp_get_max_threads()) {
     if (threads < 1) threads = 1;
     omp_set_num_threads(threads);
  }

  if(correctness) {
    arr = gen_int_array(value, n);
    erg = single_sum(arr, n);
    free(arr);
    value = 1;
  }

 // fprintf(stdout, "--- Stats ---\n\n");
 // fprintf(stdout, "n: %f, n*log(n): %f, n*log(n)^2: %f, n*n: %f\n\n",
 //	  (double) n, n*log2(n), n*log2(n)*log2(n), (double) n*n);

  if(reduction) {
    if(single)
       testReduce("Single Reduction", n, value, erg, correctness, &single_sum);
    if(all_threads) {
      for(i=1; i<=threads; i++) {
         testParaReduce("Parallel Reduction", n, value, erg, correctness, i, &reduce_scan);
      }
    } else
      testParaReduce("Parallel Reduction", n, value, erg, correctness, threads, &reduce_scan);
    fprintf(stdout, "\n");
  }

  if(recursive) {
    if(work) {
      if(single)
         testPrefixsum("Single Recursion", n, value, correctness, &single_recursive_scan);
      if(all_threads) {
        for(i=1; i<=threads; i++) {
          testParaPrefixsum("Parallel Recursion", n, value, correctness, i, &recursive_scan);
        }
      } else
         testParaPrefixsum("Parallel Recursion", n, value, correctness, threads, &recursive_scan);
      fprintf(stdout, "\n");
      if(only_sum) {
         testParaSum("Parallel Recursion Sum", n, value, erg, correctness, threads, &recursive_scan_sum);
         fprintf(stdout, "\n");
      }
    }
    if(performance) {
      if(single)
         testPrefixsum_time("Single Recursion", n, value, correctness, &single_recursive_scan_time);
      if(all_threads) {
        for(i=1; i<=threads; i++) {
          testParaPrefixsum_time("Parallel Recursion", n, value, correctness, i, &recursive_scan_time);
        }
      } else
         testParaPrefixsum_time("Parallel Recursion", n, value, correctness, threads, &recursive_scan_time);
      fprintf(stdout, "\n");
      if(only_sum) {
         testParaSum("Parallel Recursion", n, value, erg, correctness, threads, &recursive_scan_sum);
         fprintf(stdout, "\n");
      }
    }
  }

  if(iterative) {
    if(work) {
      if(single)
         testPrefixsum("Single Iterative", n, value, correctness, &single_iterative_scan);
      if(all_threads) {
        for(i=1; i<=threads; i++) {
         testParaPrefixsum("Parallel Hillis-Steele", n, value, correctness, i, &hillis_steele_scan);
        }
      } else
         testParaPrefixsum("Parallel Iterative", n, value, correctness, threads, &iterative_scan);
      fprintf(stdout, "\n");
      if(only_sum) {
         testParaSum("Parallel Iterative", n, value, erg, correctness, threads, &iterative_scan_sum);
         fprintf(stdout, "\n");
      }
    }
    if(performance) {
      if(single)
         testPrefixsum_time("Single Iterative", n, value, correctness, &single_iterative_scan_time);
      if(all_threads) {
        for(i=1; i<=threads; i++) {
         testParaPrefixsum_time("Parallel Hillis-Steele", n, value, correctness, i, &hillis_steele_scan_time);
        }
      } else
         testParaPrefixsum_time("Parallel Iterative", n, value, correctness, threads, &iterative_scan_time);
      fprintf(stdout, "\n");
      if(only_sum) {
         testParaSum("Parallel Iterative", n, value, erg, correctness, threads, &iterative_scan_sum);
         fprintf(stdout, "\n");
      }
    }
  }

  if(hills) {
    if(work) {
      if(single)
         testPrefixsum("Single Hillis-Steele", n, value, correctness, &single_hillis_steele_scan);
      if(all_threads) {
        for(i=1; i<=threads; i++) {
         testParaPrefixsum("Parallel Hillis-Steele", n, value, correctness, i, &hillis_steele_scan);
        }
      } else
         testParaPrefixsum("Parallel Hillis-Steele", n, value, correctness, threads, &hillis_steele_scan);
      fprintf(stdout, "\n");
    }
    if(performance) {
      if(single)
         testPrefixsum_time("Single Hillis-Steele", n, value, correctness, &single_hillis_steele_scan_time);
      if(all_threads) {
        for(i=1; i<=threads; i++) {
         testParaPrefixsum_time("Parallel Hillis-Steele", n, value, correctness, i, &hillis_steele_scan_time);
        }
      } else
         testParaPrefixsum_time("Parallel Hillis-Steele", n, value, correctness, threads, &hillis_steele_scan_time);
      fprintf(stdout, "\n");
    }
  }
  return 0;
}

void testReduce(char* name, STYPE n, ATYPE value, ATYPE erg, int correct, ATYPE (*sum)(ATYPE*, STYPE)) {
  ATYPE* arr;
  double time_begin, time_end;
  ATYPE e;

  arr = gen_int_array(value, n);

  time_begin = omp_get_wtime();
  e = sum(arr, n);
  time_end = omp_get_wtime();

  if(!correct || e == erg) {
    fprintf(stdout, "%lf;",(time_end - time_begin));
  } else {
    fprintf(stdout, "false %s\n", name);
  }

  free(arr);
}

void testParaReduce(char* name, STYPE n, ATYPE value, ATYPE erg, int correct, int threads, ATYPE (*sum)(ATYPE*, STYPE, int)) {
  ATYPE* arr;
  double time_begin, time_end;
  ATYPE e;

  arr = gen_int_array(value, n);

  time_begin = omp_get_wtime();
  e =  sum(arr, n, threads);
  time_end = omp_get_wtime();

  if(!correct || e == erg) {
    fprintf(stdout, "%lf;",(time_end - time_begin));
  } else {
    fprintf(stdout, "false %s\n", name);
  }
  free(arr);
}

void testPrefixsum_time (char* name, STYPE n, ATYPE value, int correct, void (*prefixsum)(ATYPE**, STYPE)) {
  ATYPE* arr;
  double time_begin, time_end;

  arr = gen_int_array(value, n);

  time_begin = omp_get_wtime();
  prefixsum(&arr, n);
  time_end = omp_get_wtime();

  if(!correct || check_correctness(arr, n)) {
    fprintf(stdout, "%lf;",(time_end - time_begin));
    fflush(stdout);
  } else {
    fprintf(stdout, "false %s\n", name);
  }

  free(arr);
}

void testParaPrefixsum_time (char* name, STYPE n, ATYPE value, int correct ,int threads, void (*prefixsum)(ATYPE**, STYPE, int)) {
  ATYPE* arr;
  double time_begin, time_end;

  arr = gen_int_array(value, n);

  time_begin = omp_get_wtime();
  prefixsum(&arr, n, threads);
  time_end = omp_get_wtime();

  if(!correct || check_correctness(arr, n)) {
    fprintf(stdout, "%lf;", (time_end - time_begin));
    fflush(stdout);
  } else {
    fprintf(stdout, "false %s\n", name);
  }

  free(arr);
}

void testPrefixsum(char* name, STYPE n, ATYPE value, int correct, void (*prefixsum)(ATYPE**, STYPE, OTYPE*)) {
  ATYPE* arr;
  double time_begin, time_end;
  OTYPE* ops = calloc(omp_get_max_threads(), sizeof(OTYPE));
  memset((void *) ops, 0, sizeof(OTYPE)*omp_get_max_threads());
  int threadnum;

  arr = gen_int_array(value, n);

  time_begin = omp_get_wtime();
  prefixsum(&arr, n, ops);
  time_end = omp_get_wtime();

  if(!correct || check_correctness(arr, n)) {
    fprintf(stdout, "%lf;",(time_end - time_begin));
    for(threadnum=0;threadnum<omp_get_max_threads();threadnum++) {
      fprintf(stdout, "%li;",ops[threadnum]);
    }
    fprintf(stdout, "\n");
  } else {
    fprintf(stdout, "false %s\n", name);
  }

  free(ops);
  free(arr);
}

void testParaPrefixsum(char* name, STYPE n, ATYPE value, int correct ,int threads, void (*prefixsum)(ATYPE**, STYPE, OTYPE*, int)) {
  ATYPE* arr;
  double time_begin, time_end;
  OTYPE* ops = calloc(omp_get_max_threads(), sizeof(OTYPE));
  memset((void *) ops, 0, sizeof(OTYPE)*omp_get_max_threads());
  int threadnum;

  arr = gen_int_array(value, n);

  time_begin = omp_get_wtime();
  prefixsum(&arr, n, ops, threads);
  time_end = omp_get_wtime();

  if(!correct || check_correctness(arr, n)) {
    fprintf(stdout, "%lf;", (time_end - time_begin));
    for(threadnum=0;threadnum<omp_get_max_threads();threadnum++) {
      fprintf(stdout, "%li;",ops[threadnum]);
    }
    fprintf(stdout, "\n");
  } else {
    fprintf(stdout, "false %s\n", name);
  }

  free(ops);
  free(arr);
}

void testParaSum(char* name, STYPE n, ATYPE value, ATYPE sum, int correct, int threads, ATYPE (*prefixsum)(ATYPE**, STYPE, OTYPE*, int)) {
  ATYPE* arr;
  ATYPE erg;
  double time_begin, time_end;
  OTYPE* ops = calloc(omp_get_max_threads(), sizeof(OTYPE));
  memset((void *) ops, 0, sizeof(OTYPE)*omp_get_max_threads());
  int threadnum;

  arr = gen_int_array(value, n);

  time_begin = omp_get_wtime();
  erg = prefixsum(&arr, n, ops, threads);
  time_end = omp_get_wtime();

  if(!correct || erg == sum) {
    fprintf(stdout, "%lf;",(time_end - time_begin));
    for(threadnum=0;threadnum<omp_get_max_threads();threadnum++)
      fprintf(stdout, "%li;",ops[threadnum]);
    fflush(stdout);
  } else {
    fprintf(stdout, "false %s: erg: %li, sum: %li -- ", name, erg, sum);
    printArray(arr, n, 0);
  }

  free(ops);
  free(arr);
}

void printArray(ATYPE * arr, ATYPE size, ATYPE i) {

  for(; i<size; i++) {
    fprintf(stdout, "%li ", arr[i]);
  }
  fprintf(stdout, "\n");
}

ATYPE* gen_int_array(ATYPE v, STYPE n) {
  ATYPE* ret = calloc(n, sizeof(ATYPE));

  for(int i = 0; i < n; i++)
     ret[i] = v;

  return ret;
}


ATYPE single_sum(ATYPE* x, STYPE n) {
  ATYPE sum = 0;

  for(STYPE i = 0; i < n;i++) {
     sum += x[i];
  }

  return sum;
}

ATYPE reduce_scan(ATYPE* x, STYPE n, int threads) {
  ATYPE sum = 0;

  #pragma omp parallel for num_threads(threads) reduction(+:sum)
  for(STYPE i = 0; i < n; i++) {
     sum += x[i];
  }

  return sum;
}

int check_correctness(ATYPE* x, STYPE n) {
  STYPE i;
  int error = 0;

  for(i=0; i < n; i++) {
    error = error || ((i+1) ^ x[i]);
  }

  return !error;
}

static void usage(void) {
   (void) fprintf(stderr, "Usage: prefixsum [-n array_length ] [ -v value ] [-t amnt_threads | -a | -m | -s ] [ -c | -p | -w ] [ -r | -i | -h | -e | -o ]\n");
   exit(EXIT_FAILURE);
}