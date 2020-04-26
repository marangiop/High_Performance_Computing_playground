// 
//The objective of this piece of work is to implement a two-dimensional lattice-based calculation that employs a two-dimensional domain decomposition and uses non-blocking communications. The developed code was applied within the context of image reconstruction from a starting image containing edges. This operation is conceptually similar to a large number of real scientific HPC calculations that solve partial differential equations using iterative algorithms such as Jacobi or Gauss-Seidel.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
#include "pgmio.h"

#define M 256 
#define N 192
#define MAXITER   1600
#define PRINTFREQ  5
#define ndims 2

double boundaryval(int i, int m);

int main (int argc, char **argv) {
  double val;
  double masterbuf[M][N];

  int i, j, iter, maxiter;
  char *filename;

  int rank, size, next, prev;

  int dims[ndims];
  dims[0] = 0;
  dims[1] = 0;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);


  //READING IN FILE
  if (rank == 0) //Rank 0 reads the input image file into a masterbuf array of size M X N with no halos
  {
    printf("Processing %d x %d image on %d processes\n", M, N, size);
    printf("Number of iterations = %d\n", MAXITER);

    filename = "edgenew256x192.pgm";

    printf("\nReading <%s>\n", filename);
    pgmread(filename, &masterbuf[0][0], M, N);
    printf("\n");
  }


  //CREATE CARTESIAN TOPOLOGY

  MPI_Comm cart_comm;

  MPI_Dims_create(size, ndims, dims);

  int period[ndims];
  int reorder;

  period[0] = 0;
  period[1] = 1; // Up and down boundaries are set to be periodic
  reorder = 0;

  MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, period, reorder, &cart_comm);
  MPI_Bcast(&masterbuf, M * N, MPI_DOUBLE, 0, cart_comm);

  double t1, t2;

    if (rank==0){
        t1 = MPI_Wtime();
    }

  int cartesian_coordinate[2];

  int MP = M / dims[0];
  int NP = N / dims[1];

  double local_array[MP][NP];
  MPI_Cart_coords(cart_comm, rank, 2, cartesian_coordinate);


  //2D DECOMPOSITION
  for (i = 0; i < MP; i++) {
    for (j = 0; j < NP; j++) {
      local_array[i][j] = masterbuf[(cartesian_coordinate[0] * MP + i)][(cartesian_coordinate[1] * NP + j)];
    }
    }
    
  double old[MP + 2][NP+ 2], new[MP + 2][NP + 2], edge[MP + 2][NP + 2];


  for (i = 1; i < MP + 1; i++) {
    for (j = 1; j < NP + 1; j++) {
      edge[i][j] = local_array[i - 1][j - 1];
    }
  }

  for (i = 0; i < MP + 2; i++) {
    for (j = 0; j < NP + 2; j++) {
      old[i][j] = 255.0;
    }
  }


  for (i = 1; i < NP + 1; i++) {
    val = boundaryval(j, N);

    old[0][j] = (int) (255.0 * (1.0 - val));
    old[M + 1][j] = (int) (255.0 * val);
  }


  int right, left, down, up;
  MPI_Cart_shift( cart_comm , 1 , 1 , &left , &right );
  MPI_Cart_shift( cart_comm , 0 , 1 , &down , &up );


  //DATATYPE FOR SENDING ROWS
  MPI_Datatype xSlice;
  MPI_Type_vector(MP, 1, NP, MPI_DOUBLE, &xSlice);

  MPI_Type_commit (& xSlice );

  MPI_Status status;
  MPI_Request request;


  //START OF ITERATIVE LOOP
  for (iter = 1; iter <= MAXITER; iter++) {


    if (iter % PRINTFREQ == 0) {
      if (rank == 0) {
        printf("Iteration %d\n", iter);
      }
    }

    for (i = 1; i < M + 1; i++) {
      old[i][0] = old[i][N];
      old[i][N + 1] = old[i][1];
    }

    //HALO SWAPS USING NON BLOCKING COMMUNICATION

    MPI_Issend(&old[1][NP+1] ,1, xSlice , up,   0, cart_comm, &request); //SENDING THE HIGHEST ARRAY ROW TO THE PROCESS THAT IS ON TOP OF ME
    MPI_Recv(&old[1][0]        ,1, xSlice , down, 0, cart_comm, &status); //WE DO NOT CONSIDER THE CORNER //MESSAGE FROM PROCESS BELOW IS RECEIVED IN THE LOWEST HALO ROW OF THE PROCESS THAT IS ON TOP
    MPI_Wait(&request,&status);

    MPI_Issend(&old[1][1]      ,1, xSlice , down, 0, cart_comm, &request );//SENDING THE LOWEST ARRAY ROW TO THE PROCESS THAT IS BELOW ME
    MPI_Recv(&old[1][NP+1]   ,1, xSlice , up,   0, cart_comm, &status ); //MESSAGE FROM PROCESS ABOVE IS RECEIVED IN THE HIGHEST HALO ROW OF THE PROCESS THAT IS BELOW //1 instead of 2 i nrecv
    MPI_Wait(&request,&status);

    MPI_Issend(&old[MP+1][1]  ,NP, MPI_DOUBLE ,right, 0, cart_comm, &request); // SENDING THE RIGHTMOST ARRAY COLUMN TO THE PROCESS THAT IS RIGHT TO ME
    MPI_Recv(&old[0][1]       ,NP, MPI_DOUBLE , left,   0, cart_comm, &status ); //MESSAGE FROM PROCESS ON THE LEFT IS RECEIVED IN THE LEFTMOST HALO COLUMN IN THE PROCESS THAT IS ON THE RIGHT
    MPI_Wait(&request,&status);

    MPI_Issend(&old[1][1]      ,NP, MPI_DOUBLE,left,  0, cart_comm, &request ); //SENDING THE LEFTMOST ARRAY COLUMN TO THE PROCESS THAT IS LEFT TO ME
    MPI_Recv(&old[MP+2][1]    ,NP, MPI_DOUBLE, right, 0, cart_comm, &status ); //MESSAGE FROM PROCESS ON THE RIGHT IS RECEIVED IN THE RIGHTMOST HALO COLUMN IN THE PROCESS THAT IS ON THE LEFT
    MPI_Wait(&request,&status);


    for (i=1;i<MP +1;i++)
    {
      for (j=1;j<NP +1;j++)
      {
        new[i][j]=0.25*(old[i-1][j]+old[i+1][j]+old[i][j-1]+old[i][j+1] - edge[i][j]);
      }
    }

    for (i=1;i<MP+1;i++)
    {
      for (j=1;j<NP+1;j++)
      {
        old[i][j]=new[i][j];
      }
    }
  }

  if (rank==0)
  {
    printf("\nFinished %d iterations\n", iter-1);
  }

  for (i=1;i<MP+1;i++)
  {
    for (j=1;j<NP+1;j++)
    {
      local_array[i-1][j-1]=old[i][j];
    }
  }

  double bigbuf[M][N];

  for (i=0; i<M; i++)
  {
    for (j=0; j<N; j++)
    {
      bigbuf[i][j]=0.0;
    }
  }

  MPI_Cart_coords(cart_comm, rank, 2, cartesian_coordinate);

  for (i=0; i<MP; i++)
  {
    for (j=0; j<NP; j++)
    {
      bigbuf[cartesian_coordinate[0]*MP+i][cartesian_coordinate[1]*NP+j]=local_array[i][j]; //Copy local picture from local_array to correct place in bigbuf
    }
  }
      
  MPI_Reduce(bigbuf, masterbuf, M*N, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); //Reduce partial images to masterbuf on rank 0 */

    if (rank==0)
    {
        t2 = MPI_Wtime();
        printf( "Elapsed time is %lf\n", t2 - t1);
    }

  //Write to file
  if (rank == 0)
  {
    filename="imagenew256x192.pgm";
    printf("\nWriting <%s>\n", filename);
    pgmwrite(filename, &masterbuf[0][0], M, N);
  }

    MPI_Finalize();
}

double boundaryval(int i, int m)
{
  double val;

  val = 2.0*((double)(i-1))/((double)(m-1));
  if (i >= m/2+1) val = 2.0-val;

  return val;
}
