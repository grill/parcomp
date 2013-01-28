#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <mpi.h>
#include <time.h>
#include <assert.h>
#include <string.h>
#include <math.h>

#define ATYPE long int

static void usage(void);
void printArray(int* arr, ATYPE size, int rank);
double radixsort(ATYPE R, ATYPE n, int size, int rank, int print_array);
double bucket_sort(ATYPE R, ATYPE n, int size, int rank, int print_array, ATYPE radix);

int main(int argc, char *argv[])
{
  int rank;
  int size;
  double time;
  double* rtime = malloc(sizeof(double));

  MPI_Init(&argc,&argv);

  // get rank and size from communicator
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

//  if (rank==0) {
//    printf("Rank %d initializing, total %d\n",rank,size);
//  }

  ATYPE R = 10; //Numbers in array
  ATYPE radix = 10;
  ATYPE n = 100; //size of array
  int r = 0;
  int print_array = 0;

  char c;

   opterr = 0;
   while ((c = getopt(argc, argv, "n:R:A:dr")) != -1 ) {
      switch (c) {
      case 'n': //array length
         if(sscanf(optarg,"%ld",&n) == EOF) {
            (void) fprintf(stderr, "%s: %s is not a valid number.\n", argv[0],
               optarg);
            usage();
         }
         break;
      case 'R': //array length
         if(sscanf(optarg,"%ld",&R) == EOF) {
            (void) fprintf(stderr, "%s: %s is not a valid number.\n", argv[0],
               optarg);
            usage();
         }
         break;
      case 'A':
         if(sscanf(optarg,"%ld",&radix) == EOF) {
            (void) fprintf(stderr, "%s: %s is not a valid number.\n", argv[0],
               optarg);
            usage();
         }
         break;
      case 'd': print_array = 1; break;
      case 'r': r = 1; break;
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

   if(r)
     time = radixsort(R, n, size, rank, print_array);
  else
     time = bucket_sort(R, n, size, rank, print_array, radix);

  MPI_Reduce(&time, rtime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  if(rank == 0) {
    fprintf(stdout, "%f", (*rtime));
  }

  free(rtime);
  MPI_Finalize();
  return 0;
}

double bucket_sort(ATYPE R, ATYPE n, int size, int rank, int print_array) {
  ATYPE localSize = n/size;
  //step 1: local bucket sort; buckets[i] = amnt of keys i
  ATYPE i;
  ATYPE k;
  ATYPE amntRecvel = 0;

  double start;
  double end;

  int* bucket = calloc(R, sizeof(int));
  int* LocB = calloc(R, sizeof(int) );
  int* AllB = calloc(R, sizeof(int));
  int* RelB = calloc(R, sizeof(int));

  int* sendelts = calloc(size, sizeof(int));
  int* sdispls = calloc(size, sizeof(int));
  int* recvelts = calloc(size, sizeof(int));
  int* rdispls = calloc(size, sizeof(int));

  memset((void* ) sendelts, 0, sizeof(int)*size);

  int* A = calloc(localSize, sizeof(int));
  int* B = calloc(localSize, sizeof(int));
  int* C;

  /* Seed the random number generator */
  srand(time(NULL));

  /* Initialize array with random numbers, from 0 to max_num */
  for(i = 0; i < localSize; i++)
    A[i] = rand() % R;
  if(print_array)
    for(i = 0; i < localSize; i++)
      fprintf(stderr,"A[%li] = %d\n",i,A[i]);

//  if(rank == 0)
//    fprintf(stdout, "Bucket Sort...\n\n");

  /* Start Timer */
  start = MPI_Wtime();

  for (i=0; i<R; i++) bucket[i] = 0;
  for (i=0; i<localSize; i++) bucket[A[i]]++;

  LocB[0] = bucket[0];
  for (i=1; i<R; i++) {
    LocB[i] = bucket[i];
    bucket[i] += bucket[i-1];
  }
  for (i=R-1; i>0; i--) bucket[i] = bucket[i-1];
  bucket[0] = 0;

  for (i=0; i<localSize; i++) B[bucket[A[i]]++] = A[i];

  free(bucket);
  //step 2: O(n + log p)
  MPI_Allreduce(LocB,AllB,R,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

  //step 3
  MPI_Exscan(LocB,RelB,R,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

  //Step 4: local exclusive prefix-sums of AllB
  for (i=1; i<R; i++) AllB[i] += AllB[i-1];
  for (i=R-1; i>0; i--) AllB[i] = AllB[i-1];
  AllB[0] = 0;

  //Now: Local element A[j] needs to go to position
  //AllB[A[j]]+RelB[A[j]]+j‘ (for A[j]>0)

  //Step 5: compute number of elements to be sent to each other process, sendelts[i], i=0,...,p-1
  for(i=0; i<localSize; i++) {
    k = (AllB[B[i]]+RelB[B[i]]++) / localSize;
    if(k != rank) {
   //   printf("%i i: %i, B[i]: %i, k: %i, localSize: %i, bigdiv: %i\n",
	//     rank, i, B[i], k, localSize, AllB[B[i]]+(RelB[B[i]])-1);
      LocB[B[i]]--;
      sendelts[k]++;
      amntRecvel++;
      if(sendelts[k] == 1)
         sdispls[k] = i;
    } else {
      A[i-amntRecvel] = B[i];
    }
  }
  free(AllB);
  free(RelB);
  C = calloc(amntRecvel, sizeof(int));

  //Step 6:
  MPI_Alltoall(sendelts,1,MPI_INT,recvelts,1,MPI_INT,MPI_COMM_WORLD);

  rdispls[0] = 0;
  for (i=1; i<size; i++) {
    rdispls[i] = recvelts[i-1] + rdispls[i-1];
  }

  //Step 7: redistribute elements
  MPI_Alltoallv(B,sendelts,sdispls,MPI_INT,C,recvelts, rdispls, MPI_INT, MPI_COMM_WORLD);

  k = 0;
  //Step 8: reorder elements from C back to A
  for(i = 0; i < amntRecvel; i++) LocB[C[i]]++;

  for (i=1; i<R; i++) LocB[i] += LocB[i-1];
  for (i=R-1; i>0; i--) LocB[i] = LocB[i-1];
  LocB[0] = 0;

  //for(i = 0; i < amntRecvel; i++)
  //  fprintf(stderr,"%i C[%d] = %d, LocB[C[%d]] = %d\n",rank,i,C[i], i, LocB[C[i]]);

  for (i=0; i<localSize-amntRecvel; i++) {
   // printf("%i LocB[A]: %i, A: %i, i: %i\n",
//	     rank, LocB[A[i]], A[i], i);
    B[LocB[A[i]]++] = A[i];
  }

 // printf("first\n");
  //printArray(B, localSize, rank);

  for (i=0; i<amntRecvel; i++) {
    //printf("%i LocB[C]: %i, C: %i, i: %i\n",
	//     rank, LocB[C[i]], C[i], i);
    B[LocB[C[i]]++] = C[i];
  }

  end = MPI_Wtime();
 
  if(print_array)
    printArray(B, localSize, rank);
//  fprintf(stdout, "%f", (end - start));
//  fprintf(stdout, "Proc %i - Time: %f\n", rank, (end - start));

  free(A);
  free(LocB);
  free(C);
  free(B);
  free(sendelts);
  free(sdispls);
  free(rdispls);
  free(recvelts);

  //Possible optimization: replace MPI_Allreduce by MPI_Bcast
  return end - start;
}

double radixsort(ATYPE R, ATYPE n, int size, int rank, int print_array, ATYPE radix) {
  int localSize = n/size;
  int i;
  int k;
  int amntRecvel;
  int exp = 1;

  double start;
  double end;

  int* bucket = calloc(radix, sizeof(int));
  int* LocB = calloc(radix, sizeof(int));
  int* AllB = calloc(radix, sizeof(int));
  int* RelB = calloc(radix, sizeof(int));

  int* sendelts = calloc(size, sizeof(int));
  int* sdispls = calloc(size, sizeof(int));
  int* recvelts = calloc(size, sizeof(int));
  int* rdispls = calloc(size, sizeof(int));

  int* A = calloc(localSize, sizeof(int));
  int* B = calloc(localSize, sizeof(int));
  int* C = calloc(localSize, sizeof(int));


  /* Seed the random number generator */
  srand(time(NULL));

  /* Initialize array with random numbers, from 0 to R */
  for(i = 0; i < localSize; i++)
    A[i] = rand() % R;
 // if(print_array)
 //   for(i = 0; i < localSize; i++)
  //    fprintf(stderr,"A[%d] = %d\n",i,A[i]);

//  if(rank == 0)
//    fprintf(stdout, "Radix Sort...\n\n");

  /* Start Timer */
  start = MPI_Wtime();

  while (R / exp > 0) {
    memset((void* ) bucket, 0, sizeof(int)*radix);
    memset((void* ) AllB, 0, sizeof(int)*radix);
    memset((void* ) RelB, 0, sizeof(int)*radix);

    //step 1: local bucket sort; buckets[i] = amnt of keys i
    for (i = 0; i < localSize; i++) bucket[(A[i] / exp) % radix]++;

    LocB[0] = bucket[0];
    for (i = 1; i < radix; i++) {
      LocB[i] = bucket[i];
      bucket[i] += bucket[i - 1];
    }
    for (i=radix-1; i>0; i--) bucket[i] = bucket[i-1];
    bucket[0] = 0;

    for (i = 0; i < localSize; i++) B[bucket[(A[i] / exp) % radix]++] = A[i];

    //step 2: O(n + log p)
    MPI_Allreduce(LocB,AllB,radix,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

    //step 3
    MPI_Exscan(LocB,RelB,radix,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

    //Step 4: local exclusive prefix-sums of AllB
    for (i=1; i<radix; i++) AllB[i] += AllB[i-1];
    for (i=radix-1; i>0; i--) AllB[i] = AllB[i-1];
    AllB[0] = 0;

    //Step 5: compute number of elements to be sent to each other process, sendelts[i], i=0,...,p-1
    amntRecvel = 0;
    memset((void* ) sendelts, 0, sizeof(int)*size);

    for(i=0; i<localSize; i++) {
      k = (AllB[(B[i] / exp) % radix]+RelB[(B[i] / exp) % radix]++) / localSize;
      if(k != rank) {
        //printf("%i i: %i, B[i]: %i, mB[i]: %i, k: %i, localSize: %i, all: %i, rel: %i\n",
	 //    rank, i, B[i], (B[i] / exp) % radix, k,
	  //     localSize, AllB[(B[i] / exp) % radix], (RelB[(B[i] / exp) % radix])-1);
        LocB[(B[i] / exp) % radix]--;
        sendelts[k]++;
        amntRecvel++;
        if(sendelts[k] == 1)
           sdispls[k] = i;
      } else {
        A[i-amntRecvel] = B[i];
      }
    }

    //Step 6:
    MPI_Alltoall(sendelts,1,MPI_INT,recvelts,1,MPI_INT,MPI_COMM_WORLD);

    rdispls[0] = 0;
    for (i=1; i<size; i++) {
      rdispls[i] = recvelts[i-1] + rdispls[i-1];
    }

    //Step 7: redistribute elements
    MPI_Alltoallv(B,sendelts,sdispls,MPI_INT,C,recvelts, rdispls, MPI_INT, MPI_COMM_WORLD);

    k = 0;
    //Step 8: reorder elements from C back to A
    for(i = 0; i < amntRecvel; i++) LocB[(C[i] / exp) % radix]++;

    for (i=1; i<radix; i++) LocB[i] += LocB[i-1];
    for (i=radix-1; i>0; i--) LocB[i] = LocB[i-1];
    LocB[0] = 0;

    //printArray(B, localSize, rank);

    for(i=0; i < rank; i++) {
      for(k=0; k < recvelts[i]; k++)
	B[LocB[(C[rdispls[i]+k] / exp) % radix]++] = C[rdispls[i]+k];
    }

    for (i=0; i<localSize-amntRecvel; i++) {
      B[LocB[(A[i] / exp) % radix]++] = A[i];
    }

    //printf("bl: %i\n", amntRecvel);
    //printArray(C, localSize, rank);
    //printArray(LocB, radix, rank);
    //printArray(B, localSize, rank);

    for(i=rank+1; i < size; i++) {
      for(k=0; k < recvelts[i]; k++)
	B[LocB[(C[rdispls[i]+k] / exp) % radix]++] = C[rdispls[i]+k];
    }
  //  for (i=0; i<amntRecvel; i++) {
  //    B[LocB[(C[i] / exp) % radix]++] = C[i];
  //  }

    for (i = 0; i < localSize; i++)
      A[i] = B[i];
    exp *= radix;
 
    if(print_array) {
      printf("\nPASS %.0f  -  ", log10(exp)/log10(radix));
      printArray(A, localSize, rank);
    }
  }

  end = MPI_Wtime();
//  fprintf(stdout, "%f", (end - start));
//  fprintf(stdout, "Proc %i - Time: %f\n", rank, (end - start));

  free(AllB);
  free(RelB);
  free(LocB);
  free(A);
  free(C);
  free(B);
  free(sendelts);
  free(sdispls);
  free(rdispls);
  free(recvelts);

  return end - start;
}

void printArray(int * arr, ATYPE size, int rank) {
  ATYPE i;

  fprintf(stdout, "Proc %i: ", rank);
  for(i=0; i<size; i++) {
    fprintf(stdout, "%d ", arr[i]);
  }
  fprintf(stdout, "\n");
}

static void usage(void) {
   (void) fprintf(stderr, "Usage: bucket-sort [-a print_array | -n size | -R nums]\n");
   exit(EXIT_FAILURE);
}
