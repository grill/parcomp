/* TUW, October 2011 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <getopt.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>


// MPI header
#include <mpi.h>

#define INT_T uint64_t
#define INT_MPI_T MPI_UINT64_T
#define HELLO 1234 // tag for control messages
#define SCAN 1337

struct opt {
    INT_T size;
    int parallel;
};

void arrayscan(INT_T A[], INT_T n, MPI_Comm comm, INT_T (*commscan)(INT_T, MPI_Comm)) ;
INT_T* gen_data(int rank, INT_T size) ;
void localscan(INT_T A[], INT_T n) ;
INT_T commscan_primitive(INT_T A, MPI_Comm comm) ;
INT_T my_commscan(INT_T A, MPI_Comm comm) ;
void parse_args(int argc, char ** argv, struct opt* args) ;
INT_T check_asc(INT_T* A, INT_T size) ;

int main(int argc, char *argv[])
{
  int rank, comm_size;
  int prev;
  char name[MPI_MAX_PROCESSOR_NAME];
  int nlen;
  INT_T size = 100000;
  double time = 0;
  struct opt args;

  MPI_Init(&argc,&argv);

  // get rank and size from communicator
  MPI_Comm_size(MPI_COMM_WORLD,&comm_size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  //printf("%d %d\n", rank, argc);
  parse_args(argc, argv, &args);
  if(args.size > 0) {
    size = args.size/comm_size;
  }

  MPI_Get_processor_name(name,&nlen);

  INT_T* A = gen_data(rank, size);
  arrayscan(A, size, MPI_COMM_WORLD, commscan_primitive);

#ifdef DEBUG
  {
    INT_T a = check_asc(A, size);
    if(a == -1) {
        printf("ascending check ok\n");
    } else {
        printf("ascending error on position %ld\n", a);
    }
    #ifdef DEBUGDEBUG
    for(INT_T i = 0; i < size; i++) {
        printf("%ld\n", A[i]);
    }
    #endif
  }
#endif
  time = - MPI_Wtime();
/*
  if (rank==0) {
    printf("Rank %d initializing, total %d\n",rank,comm_size);
  } else {
    MPI_Recv(&prev,1,MPI_INT,rank-1,HELLO,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    printf("Rank %d on %s received from %d, passing on\n",rank,name,prev);
  }
  if (rank+1<size) MPI_Send(&rank,1,MPI_INT,rank+1,HELLO,MPI_COMM_WORLD);
*/
  
  if(args.parallel)
    arrayscan(A, size, MPI_COMM_WORLD, my_commscan);
  else
    localscan(A, size);
    
  time += MPI_Wtime();

  printf("Rank %3d min: %20ld sum: %20ld time: %2lf\n", rank, A[0], A[size-1], time);
  if(rank == 0) {
    printf("Sum should be %ld\n", size*comm_size*(size*comm_size+1)/2);
  }
#ifdef DEBUGDEBUG
    for(INT_T i = 0; i < size; i++) {
        printf("%ld\n", A[i]);
    }
#endif
      
  MPI_Finalize();
  return 0;
}

INT_T* gen_data(int rank, INT_T size) {
    INT_T* ret = calloc(size, sizeof(INT_T));

    if(ret == NULL) {
        printf("gen_data: %s", strerror(errno));
        exit(1);
    }

    for(INT_T i = 0; i < size; i++) {
        ret[i] = 1;
    }

    return ret;
}

void arrayscan(INT_T A[], INT_T n, MPI_Comm comm, INT_T (*commscan)(INT_T, MPI_Comm)) {
    localscan(A, n);
    INT_T B = commscan(A[n-1], comm);
    for(INT_T i = 0; i < n; i++) {
        A[i] += B;
    }
}

INT_T my_commscan(INT_T A, MPI_Comm comm) {
    int rank, size;
    INT_T send = A;
    INT_T recv = 0;

    MPI_Comm_size(comm,&size);
    MPI_Comm_rank(comm,&rank);

    for(int k = 1; k < size; k <<= 1) {
        //receive from rank - k, send to rank + k
        //printf("sending %ld from %d to %d\n", send, rank, rank + k);
        if(rank | k) {
            //send first
            if(rank + k < size) {
                MPI_Send(&send, 1, INT_MPI_T, rank + k, SCAN, comm);
            }
            if(rank - k >= 0) {
                MPI_Recv(&recv, 1, INT_MPI_T, rank - k, SCAN, comm, MPI_STATUS_IGNORE);
            }
        } else {
            //receive first
            if(rank - k >= 0) {
                MPI_Recv(&recv, 1, INT_MPI_T, rank - k, SCAN, comm, MPI_STATUS_IGNORE);
            }
            if(rank + k < size) {
                MPI_Send(&send, 1, INT_MPI_T, rank + k, SCAN, comm);
            }
        }
        send += recv;
        recv = 0;
    }
    return send - A;
}

INT_T commscan_primitive(INT_T A, MPI_Comm comm) {
    int rank, size;
    INT_T prev = 0;


    MPI_Comm_size(comm,&size);
    MPI_Comm_rank(comm,&rank);

    if(rank!=0) {
        MPI_Recv(&prev,1, INT_MPI_T, rank-1, SCAN, comm, MPI_STATUS_IGNORE);
        A += prev;
    }
    if (rank+1<size) MPI_Send(&A,1, INT_MPI_T, rank+1, SCAN, comm);

    return prev;
}

void localscan(INT_T A[], INT_T n) {
    for(INT_T i = 1; i < n; i++) {
        A[i] += A[i-1];
    }
}

void parse_args(int argc, char ** argv, struct opt* args) {
    args->size = -1;
    args->parallel = 1;
    
    char c;

    while ((c = getopt (argc, argv, "ls:")) != -1)
     switch (c)
       {
       case 's':
         args->size = strtol(optarg, NULL, 16);
         break;
       case 'l':
         args->parallel = 0;
         break;
       case '?':
         if (optopt == 's')
           fprintf (stderr, "Option -%c requires an argument.\n", optopt);
         else if (isprint (optopt))
           fprintf (stderr, "Unknown option `-%c'.\n", optopt);
         else
           fprintf (stderr,
                    "Unknown option character `\\x%x'.\n",
                    optopt);
         return;
       default:
         abort ();
       }
}

INT_T check_asc(INT_T* A, INT_T size) {
    INT_T start = A[0];

    for(INT_T i = 0; i < size; i++) {
        if(A[i] != start + i)
            return i;
    }
    return -1;
}
