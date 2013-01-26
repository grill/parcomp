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
#define STENCIL 4242
#define X(A,i,j) A[(i)+(j)*(n+2)]

struct opt {
    INT_T m;
    INT_T n;
    INT_T r;
    INT_T c;
    int iterations;
    int parallel;
};

void parse_args(int argc, char ** argv, struct opt* args) ;
void gen_data(INT_T *A, INT_T n, INT_T m) ;
INT_T* compute(INT_T *A, INT_T *B, int iterations, INT_T n, INT_T m, INT_T r, INT_T c, MPI_Comm comm, MPI_Datatype column_t) ;
void synchronize(INT_T* A, INT_T n, INT_T m, INT_T r, INT_T c, MPI_Comm comm, MPI_Datatype column_t) ;
void gen_data(INT_T *A, INT_T n, INT_T m) ;
int check(INT_T *A, INT_T n, INT_T m, int iterations) ;
int rank2rid(int rank, int r, int c);
int rank2cid(int rank, int r, int c);
int getrank(int rid, int cid, int r, int c);

int rank2rid(int rank, int r, int c) {
    return rank % r;
}
int rank2cid(int rank, int r, int c) {
    return rank / r;
}
int getrank(int rid, int cid, int r, int c) {
    if(cid < 0) {
        cid += c;
    }
    cid %= c;
    if(rid < 0) {
        rid += r;
    }
    rid %= r;

    return (cid * r) + rid;
}

int main(int argc, char ** argv) {
  int comm_size, rank;
  INT_T* A,* B, *C;
  struct opt args;
  double time;
  INT_T n,m;
  MPI_Datatype column_t;

  MPI_Init(&argc,&argv);

  // get rank and size from communicator
  MPI_Comm_size(MPI_COMM_WORLD,&comm_size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  parse_args(argc, argv, &args);

  if(args.c * args.r != comm_size) {
    printf("only have %d processes, should have (%ld*%ld=) %ld\n", comm_size, args.r, args.c, args.r*args.c);
    MPI_Finalize();
    return 1;
  } else if ((args.r != 1 && args.r % 2 != 0) || (args.c != 1 && args.c % 2 != 0)) {
    printf("r and c must be 1 or even\n");
    MPI_Finalize();
    return 1;
  } else if (args.n % args.r != 0 || args.m % args.c != 0) {
    printf("n must be divisable by r and m must be divisable by c\n");
    MPI_Finalize();
    return 1;
  }
  
  if(rank == 0) {
    printf("stencil computation, n=%ld, m=%ld, r=%ld, c=%ld, i=%d\n", args.n, args.m, args.r, args.c, args.iterations);
  }

  n = args.n/args.r;
  m = args.m/args.c;

  A = calloc((n+2)*(m+2),sizeof(INT_T));
  B = calloc((n+2)*(m+2),sizeof(INT_T));

  gen_data(A, n, m);

  MPI_Type_vector(m,1,n+2,INT_MPI_T,&column_t);
  MPI_Type_commit(&column_t);

  MPI_Barrier(MPI_COMM_WORLD);
  time = - MPI_Wtime();

  C = compute(A, B, args.iterations, n, m, args.r, args.c, MPI_COMM_WORLD, column_t);

  time += MPI_Wtime();

  if(check(C, n, m, args.iterations)) {
    //printf("check successfull!\n");
    printf("%d: %lf\n", rank, time);
  } else {
    //printf("check failed!\n");
    printf("F: %d: %lf\n", rank, time);
  }

#ifdef DEBUGDEBUG
  synchronize(A,n,m,1,1,MPI_COMM_WORLD);
  for(INT_T i = 0; i <= n+1; i++) {
    for(INT_T j = 0; j <= m+1; j++) {
      printf("%ld ", X(C,i,j));
    }
    printf("\n");
  }
#endif

  MPI_Finalize();
  return 0;
}

INT_T* compute(INT_T *A, INT_T *B, int iterations, INT_T n, INT_T m, INT_T r, INT_T c, MPI_Comm comm, MPI_Datatype column_t) {
    INT_T *h;

    while(iterations > 0) {
        synchronize(A, n, m, r, c, comm, column_t);
        
        for(INT_T i = 1; i <= n; i++) {
            for(INT_T j = 1; j <= m; j++) {
                //if(i == 1 && j == 1)
                //    printf("computing ( %ld + %ld + %ld + %ld ) / 4\n", X(A,i-1,j),
                X(B,i,j) = (X(A,i-1,j) + X(A,i,j-1) + X(A,i+1,j) + X(A,i,j+1)) / 4;
            }
        }

        h = A;
        A = B;
        B = h;

        iterations--;
    }

    return A;
}


void synchronize(INT_T* A, INT_T n, INT_T m, INT_T r, INT_T c, MPI_Comm comm, MPI_Datatype column_t) {
    int rank, other;
    int rid, cid;
    
    MPI_Comm_rank(comm, &rank);
    rid = rank2rid(rank, r, c);
    cid = rank2cid(rank, r, c);

    //rows of processors must synchronize up and down
    if(r == 1) {
        memcpy(A + (n+2)*(m+1)+1, A + n + 3, n*sizeof(INT_T));
        memcpy(A + 1, A + (n+2)*(m)+1, n*sizeof(INT_T));
    } else {
        if(rid % 2 == 0) {
            //up
            other = getrank(rid - 1, cid, r, c);
            MPI_Send(A + n + 3, n, INT_MPI_T, other, STENCIL, comm);
            MPI_Recv(A + 1, n, INT_MPI_T, other, STENCIL, comm, MPI_STATUS_IGNORE);
            //down
            other = getrank(rid + 1, cid, r, c);
            MPI_Recv(A + (n+2)*(m+1)+1, n, INT_MPI_T, other, STENCIL, comm, MPI_STATUS_IGNORE);
            MPI_Send(A + (n+2)*m+1, n, INT_MPI_T, other, STENCIL, comm);
        } else {
            //down
            other = getrank(rid + 1, cid, r, c);
            MPI_Recv(A + (n+2)*(m+1)+1, n, INT_MPI_T, other, STENCIL, comm, MPI_STATUS_IGNORE);
            MPI_Send(A + (n+2)*m+1, n, INT_MPI_T, other, STENCIL, comm);
            //up
            other = getrank(rid - 1, cid, r, c);
            MPI_Send(A + n + 3, n, INT_MPI_T, other, STENCIL, comm);
            MPI_Recv(A + 1, n, INT_MPI_T, other, STENCIL, comm, MPI_STATUS_IGNORE);
        }
    }
    //columns of processors must synchronize left and right
    if(c == 1) {
        for(INT_T i = 1; i <= m; i++) {
            X(A,0,i) = X(A,n,i);
            X(A,n+1,i) = X(A,1,i);
        }
    } else {
        if(cid % 2 == 0) {
            //left
            other = getrank(rid, cid - 1, r, c);
            MPI_Send(A + n + 3, 1, column_t, other, STENCIL, comm);
            MPI_Recv(A + n + 2, 1, column_t, other, STENCIL, comm, MPI_STATUS_IGNORE);
            //right
            other = getrank(rid, cid + 1, r, c);
            MPI_Recv(A + 2*n + 3, 1, column_t, other, STENCIL, comm, MPI_STATUS_IGNORE);
            MPI_Send(A + 2*n + 2, 1, column_t, other, STENCIL, comm);
        } else {
            //right
            other = getrank(rid, cid + 1, r, c);
            MPI_Recv(A + 2*n + 3, 1, column_t, other, STENCIL, comm, MPI_STATUS_IGNORE);
            MPI_Send(A + 2*n + 2, 1, column_t, other, STENCIL, comm);
            //left
            other = getrank(rid, cid - 1, r, c);
            MPI_Send(A + n + 3, 1, column_t, other, STENCIL, comm);
            MPI_Recv(A + n + 2, 1, column_t, other, STENCIL, comm, MPI_STATUS_IGNORE);
        }
    }
}
    
/*
void synchronize(INT_T* A, INT_T n, INT_T m, INT_T r, INT_T c, MPI_Comm comm, MPI_Datatype column_t) {
    if(r == 1 && c == 1) {
        memcpy(A + (n+1)*(m+2), A + m + 3, m*sizeof(INT_T));
        memcpy(A + 1, A + (n)*(m+2), m*sizeof(INT_T));
        for(INT_T i = 1; i <= n; i++) {
            X(A,i,0) = X(A,i,m);
            X(A,i,m+1) = X(A,i,1);
        }
    }
    else {
        int rank, other;

        MPI_Comm_rank(comm,&rank);

        //r-1: rank - 1
        //r+1: rank + 1
        if(r > 1) {//could be that rows/columns are switched
            if(rank % r == 0) {
                other = rank + r - 1;
            } else {
                other = rank - 1;
            }
            if(rank % 2 == 0) {
                MPI_Send(A+m+3, 1, column_t, other, STENCIL, comm);
                MPI_Recv(A+m+2, 1, column_t, other, STENCIL, comm, MPI_STATUS_IGNORE);
            } else {
                MPI_Recv(A+m+2, 1, column_t, other, STENCIL, comm, MPI_STATUS_IGNORE);
                MPI_Send(A+m+3, 1, column_t, other, STENCIL, comm);
            }
            if(rank % r == (r-1)) {
                other = rank - r + 1;
            } else {
                other = rank + 1;
            }
            if(rank % 2 == 0) {
                MPI_Send(A+2*m+3, 1, column_t, other, STENCIL, comm);
                MPI_Recv(A+2*m+4, 1, column_t, other, STENCIL, comm, MPI_STATUS_IGNORE);
            } else {
                MPI_Recv(A+2*m+4, 1, column_t, other, STENCIL, comm, MPI_STATUS_IGNORE);
                MPI_Send(A+2*m+3, 1, column_t, other, STENCIL, comm);
            }
        } else {
            for(INT_T i = 1; i <= n; i++) {
                X(A,i,0) = X(A,i,m);
                X(A,i,m+1) = X(A,i,1);
            }
        }
        //c-1: rank - r
        //c+1: rank + r
        if(c > 1) {
            if(rank / r == 0) {//something is wrong here
                other = (rank + c*r - 1) % rank;
            } else {
                other = rank - r;
            }
            if(rank % 2 == 0) {//and here
                MPI_Send(A+m+3, m, INT_MPI_T, other, STENCIL, comm);
                MPI_Recv(A+1, m, INT_MPI_T, other, STENCIL, comm, MPI_STATUS_IGNORE);
            } else {
                MPI_Recv(A+1, m, INT_MPI_T, other, STENCIL, comm, MPI_STATUS_IGNORE);
                MPI_Send(A+m+3, m, INT_MPI_T, other, STENCIL, comm);
            }
            if(rank % c == (c-1)) {//here too
                other = rank - c + 1;
            } else {
                other = rank + 1;
            }
            if(rank % 2 == 0) {//aaaaaand here as well
                MPI_Send(A+(n+1)*(m+2), m, INT_MPI_T, other, STENCIL, comm);
                MPI_Recv(A+n*(m+2), m, INT_MPI_T, other, STENCIL, comm, MPI_STATUS_IGNORE);
            } else {
                MPI_Recv(A+n*(m+2), m, INT_MPI_T, other, STENCIL, comm, MPI_STATUS_IGNORE);
                MPI_Send(A+(n+1)*(m+2), m, INT_MPI_T, other, STENCIL, comm);
            }
        }
    }
}
*/

void gen_data(INT_T *A, INT_T n, INT_T m) {
    for(INT_T i = 1; i <= n; i++) {
        for(INT_T j = 1; j <= m; j++) {
           X(A,i,j) = 128*((i+j) % 2);
        }
    }
}

int check(INT_T *A, INT_T n, INT_T m, int iterations) {
    for(INT_T i = 1; i <= n; i++) {
        for(INT_T j = 1; j <= m; j++) {
            if(X(A,i,j) != ((i+j+iterations) % 2)*128) {
                printf("check failed at %ld, %ld; is %ld, should be %ld\n", i, j, X(A,i,j), (i+j+iterations)%2);
                return 0;
            }
        }
    }
    return 1;
}

void parse_args(int argc, char ** argv, struct opt* args) {
    args->n = 1024;
    args->m = 1024;
    args->r = 1;
    args->c = 1;
    args->iterations = 1;
    
    char c;

    while ((c = getopt (argc, argv, "n:m:c:r:i:")) != -1)
     switch (c)
       {
       case 'n':
         args->n = strtol(optarg, NULL, 16);
         break;
       case 'm':
         args->m = strtol(optarg, NULL, 16);
         break;
       case 'c':
         args->c = strtol(optarg, NULL, 16);
         break;
       case 'r':
         args->r = strtol(optarg, NULL, 16);
         break;
       case 'i':
         args->iterations = strtol(optarg, NULL, 10);
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
