// Daxypy in parallel
#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

#define PDAXPY 0

// Overwrites vector v with scalar1*v + scalar2*w
void daxpy(double *v, const double *w, double scalar1, double scalar2, int n)
{

  /* To parallelize - broadcast the scalars scalar1 and scalar2,
   * then split up the loop computation for the portions of v and w 
   * that you need, and combine the result at the end.
   */
   int i;

#if PDAXPY == 1
   int rank, nprocs, numrows;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

   numrows = n/nprocs;
   double temp[numrows];

   for (i = 0; i < numrows; ++i)
   {
     temp[i] = v[i*numrows + rank]*scalar1 + w[i*numrows + rank]*scalar2;
   }

   MPI_Gather(temp, numrows, MPI_DOUBLE, v, numRows, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#else
  
  for(i = 0; i < n; i++)
  {
    v[i] = v[i]*scalar1 + w[i]*scalar2;
  }
}


int main( int argc, char* argv[] )
{
	double vec1[4] = {1.0, 2.0, 3.0, 4.0};
	double vec2[4] = {1.0, 1.0, 1.0, 1.0};
	int size = 4;
	double a = 1.0;
	double b = -2.0;
  int rank, nprocs;

	MPI_Init( &argc, &argv );

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	daxpy(&vec1, &vec2, a, b, size);

  if(rank == 0){}
  	printf("Result of daxpy is:\n{ ");
  	for (int i = 0; i < size; ++i)
  	{
  		printf("%f\n", v[i]);
  	}
  	printf("}\n");
  }

	MPI_Finalize();

	return 0;
}