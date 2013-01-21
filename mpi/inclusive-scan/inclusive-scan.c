/* TUW, October 2011 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>


// MPI header
#include <mpi.h>

#define INT_T uint64_t
#define INT_MPI_T MPI_UINT64_T
#define HELLO 1234 // tag for control messages
#define SCAN 1337

void arrayscan(INT_T A[], INT_T n, MPI_Comm comm, INT_T (*commscan)(INT_T, MPI_Comm)) ;
INT_T* gen_data(int rank, INT_T size) ;
void localscan(INT_T A[], INT_T n) ;
INT_T commscan_primitive(INT_T A, MPI_Comm comm) ;
INT_T my_commscan(INT_T A, MPI_Comm comm) ;

int main(int argc, char *argv[])
{
  int rank, comm_size;
  int prev;
  char name[MPI_MAX_PROCESSOR_NAME];
  int nlen;
  INT_T size = 100000;
  double time = 0;

  MPI_Init(&argc,&argv);

  // get rank and size from communicator
  MPI_Comm_size(MPI_COMM_WORLD,&comm_size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  //printf("%d %d\n", rank, argc);

  MPI_Get_processor_name(name,&nlen);

  INT_T* A = gen_data(rank, size);
  arrayscan(A, size, MPI_COMM_WORLD, commscan_primitive);
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
  
  arrayscan(A, size, MPI_COMM_WORLD, my_commscan);
  time += MPI_Wtime();

  printf("Rank %3d sum: %20ld time: %2lf\n", rank, A[size-1], time);
      
  MPI_Finalize();
  return 0;
}

INT_T* gen_data(int rank, INT_T size) {
    INT_T* ret = calloc(size, sizeof(INT_T));

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
