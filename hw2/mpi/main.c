/* UCSB CS240A, Winter Quarter 2014
 * Main and supporting functions for the Conjugate Gradient Solver on a 5-point stencil
 *
 * NAMES:
 * PERMS:
 * DATE:
 */
#include "mpi.h"
#include "hw2harness.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define DEBUG 1 
#define PDDOT 1
#define PDAXPY 1
#define PMATVEC 1

double* load_vec( char* filename, int* k );
void save_vec( int k, double* x );
double *cgsolve(double *x, int* i, double* norm, int n);
double ddot(double *v, double *w, int n);
void matvec(double *v, double *w, int n, int niters);
void daxpy(double *v, double *w, double alpha, double beta, int n);

int rank, nprocs;

int main( int argc, char* argv[] ) {
	int writeOutX = 0;
	int n, k;
	int maxiterations = 1000;
	int niters = 0;
	double norm;
	double* b;
	double* x;
	double time;
	double t1, t2;
	int correct;

	MPI_Init( &argc, &argv );
	
	// Read command line args.
	// 1st case runs model problem, 2nd Case allows you to specify your own b vector
	if ( argc == 3 ) {
		k = atoi( argv[1] );
		n = k*k;
		// each processor calls cs240_getB to build its own part of the b vector!
	} else if  ( !strcmp( argv[1], "-i" ) && argc == 4 ) {
		b = load_vec( argv[2], &k );
	} else {
		printf( "\nCGSOLVE Usage: \n\t"
			"Model Problem:\tmpirun -np [number_procs] cgsolve [k] [output_1=y_0=n]\n\t"
			"Custom Input:\tmpirun -np [number_procs] cgsolve -i [input_filename] [output_1=y_0=n]\n\n");
		exit(0);
	}
	writeOutX = atoi( argv[argc-1] ); // Write X to file if true, do not write if unspecified.

	// Start Timer
	t1 = MPI_Wtime();
	
	// CG Solve here!
	double x_initial[n];
	x = x_initial;
	cgsolve(x, &niters, &norm, n);
	
	// End Timer
	t2 = MPI_Wtime();
	
	if ( writeOutX ) {
		save_vec( k, x );
	}
		
	// Output
	if(rank == 0)
  {
  	printf( "Problem size (k): %d\n",k);
  	if(niters > 0) {
  	  printf( "Norm of the residual after %d iterations: %lf\n", niters, norm);
  	}
  	printf( "Elapsed time during CGSOLVE: %lf\n", t2-t1);

  	correct = cs240_verify(x, k, t2-t1);

  	printf("Correct: %d\n", correct);
  }
  else
  {
    printf("-only processor 0 prints results\n");
  }
	
	// Deallocate 
	
  //if(niters > 0)
  //  free(b);

	//if(niters > 0)
	//  free(x);
	
	MPI_Finalize();
	
	return 0;
}

/*
 * Team functions
 *
 */

double *cgsolve(double *x, int *iter, double *norm, int n)
{
  int niters = 0;   // number of iterations
  int possmax = 5*sqrt(n);  // possible max. iterations
  int MAXITERS = (1000 > possmax) ? 1000 : possmax;
  double TARGRES = 1.0e-6;  // target residual
  double relres;  // relative residual
  //double *x;  // vector that we are solving for
  double r[n];
  double d[n];  // direction
  double Ad[n];
  double alpha, beta;
  double rtr;
  double rtrold;
  double normb;
  int i;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  for(i = 0; i < n; ++i)
  {
    x[i] = 0;
    r[i] = d[i] = cs240_getB(i, n);
  }
 
  normb = sqrt(ddot(r, r, n));
  rtr = ddot(r, r, n);
  relres = 1.0;

  while(relres > TARGRES && niters < MAXITERS)
  {
if(DEBUG) printf("-proc %d in cgsolve loop %d\n", rank, niters);
    ++niters;
    matvec(Ad, d, n, niters);
    alpha = rtr / ddot(d, Ad, n);
    daxpy(x, d, 1, alpha, n);
    daxpy(r, Ad, 1, -alpha, n);
    rtrold = rtr;
    rtr = ddot(r, r, n);
    beta = rtr / rtrold;
    daxpy(d, r, beta, 1, n);
    relres = sqrt(rtr) / normb;
    MPI_Bcast(&relres, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }

  printf("niters is %d on proc %d\n", niters, rank);

  // returning values to be printed after a run
  *iter = niters;
  *norm = relres;

  return x;
}

double ddot(double *v, double *w, int n)
{
  double prod = 0.0;
  int i;
#if PDDOT == 1
  double prodsum = 0.0;
  int cellsperproc = n/nprocs;
  double subset_v[cellsperproc];
  double subset_w[cellsperproc];

  MPI_Scatter(v, cellsperproc, MPI_DOUBLE, subset_v, cellsperproc, MPI_DOUBLE,
    0, MPI_COMM_WORLD);
  MPI_Scatter(w, cellsperproc, MPI_DOUBLE, subset_w, cellsperproc, MPI_DOUBLE,
    0, MPI_COMM_WORLD);

  for(i = 0; i < cellsperproc; ++i)
  {
    prod += subset_v[i] * subset_w[i];
  }

  MPI_Reduce(&prod, &prodsum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  return prodsum;
#else
  for(i = 0; i < n; ++i)
  {
    prod += v[i] * w[i];
  }

  return prod;
#endif
}

// Overwrites vector v with scalar1*v + scalar2*w
void daxpy(double *v, double *w, double scalar1, double scalar2, int n)
{
  int i;
#if PDAXPY == 1
  double localscalar1, localscalar2;
  int cellsperproc = n/nprocs;
  double subset_v[cellsperproc];
  double subset_w[cellsperproc];

  MPI_Scatter(v, cellsperproc, MPI_DOUBLE, subset_v, cellsperproc, MPI_DOUBLE,
    0, MPI_COMM_WORLD);
  MPI_Scatter(w, cellsperproc, MPI_DOUBLE, subset_w, cellsperproc, MPI_DOUBLE,
    0, MPI_COMM_WORLD);
  if(rank == 0)
  {
    localscalar1 = scalar1;
    localscalar2 = scalar2;
  }
  MPI_Bcast(&localscalar1, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&localscalar2, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  for(i = 0; i < cellsperproc; ++i)
  {
    subset_v[i] = subset_v[i]*localscalar1 + subset_w[i]*localscalar2;
  }

  if(rank == 0)
  {
    MPI_Gather(subset_v, cellsperproc, MPI_DOUBLE, v, cellsperproc, MPI_DOUBLE,
      0, MPI_COMM_WORLD);
  }
  else
  {
    MPI_Gather(subset_v, cellsperproc, MPI_DOUBLE, NULL, cellsperproc,
      MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }
#else
  for(i = 0; i < n; i++)
  {
    v[i] = v[i]*scalar1 + w[i]*scalar2;
  }
#endif
}

void matvec(double *v, double *w, int n, int niters)
{
  int k = (int)sqrt(n);
  int i;
  int r, s;   //row, column
#if PMATVEC == 1
  int cellsperproc = n / nprocs;
  int rowsperproc = cellsperproc / k;
  double subset_w[cellsperproc];
  double subset_v[cellsperproc];
  int upneighbor, downneighbor;
  double neighborcell;
  MPI_Status status;

  MPI_Scatter(w, cellsperproc, MPI_DOUBLE, subset_w, cellsperproc, MPI_DOUBLE,
    0, MPI_COMM_WORLD);

  upneighbor = (rank - 1 >= 0) ? rank - 1 : nprocs - 1;
  downneighbor = (rank + 1) % nprocs;

  for(i = 0; i < cellsperproc; ++i)
  {
    subset_v[i] = 0.0;
  }

  // In order to synchronize send/recvs, even ranks will work from bottom up
  // while odd ranks will work from top down. This way, neighbors will be
  // exchanging cells in unison.
  if(rank % 2 == 0) // even procs
  {
    for(r = rowsperproc - 1; r >= 0; --r)
    {
      for(s = 0; s < k; ++s)
      {
        i = (r * k) + s;

        subset_v[i] = 4 * subset_w[i];

        if(rank != 0 || r != 0)
        {
if(DEBUG) printf("-in iter %d proc %d talking to proc %d\n", niters, rank, upneighbor);
          MPI_Send(&subset_w[i], 1, MPI_DOUBLE, upneighbor, 1, MPI_COMM_WORLD);
          MPI_Recv(&neighborcell, 1, MPI_DOUBLE, upneighbor, 2,
            MPI_COMM_WORLD, &status);
          subset_v[i] -= neighborcell;
        }
        if(s != 0)
          subset_v[i] -= subset_w[i-1];
        if(s != k-1)
          subset_v[i] -= subset_w[i+1];
        if(rank != nprocs - 1 || r != rowsperproc-1)
        {
if(DEBUG) printf("-in iter %d proc %d talking to proc %d\n", niters, rank, downneighbor);
          MPI_Send(&subset_w[i], 1, MPI_DOUBLE, downneighbor, 1, 
            MPI_COMM_WORLD);
          MPI_Recv(&neighborcell, 1, MPI_DOUBLE, downneighbor, 2,
            MPI_COMM_WORLD, &status);
          subset_v[i] -= neighborcell;
        }
      }
    }
  }
  else    // odd procs
  {
    for(r = 0; r < rowsperproc; ++r)
    {
      for(s = 0; s < k; ++s)
      {
        i = (r * k) + s;

        subset_v[i] = 4 * subset_w[i];

        if(rank != 0 || r != 0)
        {
if(DEBUG) printf("-in iter %d proc %d talking to proc %d\n", niters, rank, upneighbor);
          MPI_Recv(&neighborcell, 1, MPI_DOUBLE, upneighbor, 1,
            MPI_COMM_WORLD, &status);
          MPI_Send(&subset_w[i], 1, MPI_DOUBLE, upneighbor, 2, MPI_COMM_WORLD);
          subset_v[i] -= neighborcell;
        }
        if(s != 0)
          subset_v[i] -= subset_w[i-1];
        if(s != k-1)
          subset_v[i] -= subset_w[i+1];
        if(rank != nprocs - 1 || r != rowsperproc-1)
        {
if(DEBUG) printf("-in iter %d proc %d talking to proc %d\n", niters, rank, downneighbor);
          MPI_Recv(&neighborcell, 1, MPI_DOUBLE, downneighbor, 1, 
            MPI_COMM_WORLD, &status);
          MPI_Send(&subset_w[i], 1, MPI_DOUBLE, downneighbor, 2,
            MPI_COMM_WORLD);
          subset_v[i] -= neighborcell;
        }
      }
    }
  }

  if(rank == 0)
  {
    MPI_Gather(subset_v, cellsperproc, MPI_DOUBLE, v, cellsperproc, 
      MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }
  else
  {
    MPI_Gather(subset_v, cellsperproc, MPI_DOUBLE, NULL, cellsperproc,
      MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }
#else
  for(i = 0; i < n; ++i)
  {
    v[i] = 0.0;
  }

  for(r = 0; r < k; ++r)
  {
    for(s = 0; s < k; ++s)
    {
      i = (r * k) + s;
      v[i] = 4 * w[i];

      if(r != 0)
        v[i] -= w[i-k];
      if(s != 0)
        v[i] -= w[i-1];
      if(s != k-1)
        v[i] -= w[i+1];
      if(r != k-1)
        v[i] -= w[i+k];
    }
  }
#endif
}


/*
 * Supporting Functions
 *
 */

// Load Function
// NOTE: does not distribute data across processors
double* load_vec( char* filename, int* k ) {
	FILE* iFile = fopen(filename, "r");
	int nScan;
	int nTotal = 0;
	int n;
	
	if ( iFile == NULL ) {
		printf("Error reading file.\n");
		exit(0);
	}
	
	nScan = fscanf( iFile, "k=%d\n", k );
	if ( nScan != 1 ) {
		printf("Error reading dimensions.\n");
		exit(0);
	}
	
	n = (*k)*(*k);
	double* vec = (double *)malloc( n * sizeof(double) );
	
	do {
		nScan = fscanf( iFile, "%lf", &vec[nTotal++] );
	} while ( nScan >= 0 );
	
	if ( nTotal != n+1 ) {
		printf("Incorrect number of values scanned n=%d, nTotal=%d.\n",n,nTotal);
		exit(0);
	}
	
	return vec;
}

// Save a vector to a file.
void save_vec( int k, double* x ) { 
	FILE* oFile;
	int i;
	oFile = fopen("xApprox.txt","w");
	
	fprintf( oFile, "k=%d\n", k );
	
	for (i = 0; i < k*k; i++) { 
    	fprintf( oFile, "%lf\n", x[i]);
 	} 

	fclose( oFile );
}

