/* UCSB CS240A, Winter Quarter 2014
 * Main and supporting functions for the Conjugate Gradient Solver on a 5-point stencil
 *
 * NAMES: Kyle Jorgensen and Jasen Hall
 * PERMS: 4165916 and 8408742
 * DATE: 28 January 2014
 */
#include "mpi.h"
#include "hw2harness.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define DEBUG 0 
#define PDDOT 1
#define PDAXPY 1
#define PMATVEC 1

#define MAXITERSLIMIT 100

double* load_vec( char* filename, int* k );
void save_vec( int k, double* x );
double *cgsolve(double *x, int* i, double* norm, int n);
double ddot(double *v, double *w, int n);
void matvec(double *v, double *w, int n);
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
  int MAXITERS = (MAXITERSLIMIT > possmax) ? MAXITERSLIMIT : possmax;
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
//if(DEBUG) printf("-proc %d in cgsolve loop %d\n", rank, niters);
    ++niters;
    matvec(Ad, d, n);
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

if(DEBUG) printf("niters is %d on proc %d\n", niters, rank);

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
  // Broadcast the scalars from proc 0 so that everyone is computing with the 
  // correct values of alpha and beta
  MPI_Bcast(&localscalar1, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&localscalar2, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  for(i = 0; i < cellsperproc; ++i)
  {
    subset_v[i] = subset_v[i]*localscalar1 + subset_w[i]*localscalar2;
  }

  // Processor 0 receives all the subset_v data from everyone, and then puts it in v
  if(rank == 0)
  {
    MPI_Gather(subset_v, cellsperproc, MPI_DOUBLE, v, cellsperproc, MPI_DOUBLE,
      0, MPI_COMM_WORLD);
  }
  else  // Everyone else just sends its subset_v data to proc 0
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

void matvec(double *v, double *w, int n)
{
  int k = (int)sqrt(n);
  int i;
  int r, s;   //row, column
#if PMATVEC == 1
  int cellsperproc = n / nprocs;
  int rowsperproc = cellsperproc / k;

  // These subset vectors store an extra 2*k amount of data because they store
  // one ghost row above and below their own set of rows.  
  double subset_v[cellsperproc + (2*k)]; 
  double subset_w[cellsperproc + (2*k)]; 
  int upneighbor, downneighbor;
  MPI_Status status;

  // Scatter the w vector on proc 0 to each processor's subset_w vector
  MPI_Scatter(w, cellsperproc, MPI_DOUBLE, subset_w + k, cellsperproc,
    MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // figure out comm partners for this processor
  // top and bottom processors will exchange data they don't use
  upneighbor = (rank - 1 >= 0) ? rank - 1 : nprocs - 1;
  downneighbor = (rank + 1) % nprocs;

  if(nprocs > 1){

    // even send/recv up, odd recv/send down
    if(rank % 2 == 0) // even
      {
      	// send first row of real data (+ k)
      	MPI_Send(subset_w + k, k, MPI_DOUBLE, upneighbor, 1, MPI_COMM_WORLD);
      	// recv first row of ghost data (+ 0)
      	MPI_Recv(subset_w, k, MPI_DOUBLE, upneighbor, 2, MPI_COMM_WORLD,
      		 &status);
      }
    else  // odd
      {
      	// recv last row of ghost data (+ cellsperproc + k)
      	MPI_Recv(subset_w + (cellsperproc + k), k, MPI_DOUBLE, downneighbor,
      		 1, MPI_COMM_WORLD, &status);
      	// send last row of real data (+ cellsperproc)
      	MPI_Send(subset_w + cellsperproc, k, MPI_DOUBLE, downneighbor, 2,
      		 MPI_COMM_WORLD);
      }
    
    // even send/recv down, odd recv/send up
    if(rank % 2 == 0) // even
      {
      	// send last row of real data (+ cellsperproc)
      	MPI_Send(subset_w + cellsperproc, k, MPI_DOUBLE, downneighbor, 3,
      		 MPI_COMM_WORLD);
      	// recv last row of ghost data (+ cellsperproc + k)
      	MPI_Recv(subset_w + (cellsperproc + k), k, MPI_DOUBLE, downneighbor,
      		 4, MPI_COMM_WORLD, &status);
      }
    else  // odd
      {
      	// recv first row of ghost data (+ 0)
      	MPI_Recv(subset_w, k, MPI_DOUBLE, upneighbor, 3, MPI_COMM_WORLD,
      		 &status);
      	// send first row of real data (+ k)
      	MPI_Send(subset_w + k, k, MPI_DOUBLE, upneighbor, 4, MPI_COMM_WORLD);
      }
  }

// if(DEBUG && (niters % 2 == 0)) {
// if(rank == 0) {
// printf("%d %d - %f\n", niters, rank, *(subset_w));
// printf("%d %d + %f\n", niters, rank, *(subset_w+k));
// } else if(rank == nprocs - 1) {
// printf("%d %d - %f\n", niters, rank, *(subset_w+cellsperproc));
// printf("%d %d + %f\n", niters, rank, *(subset_w+cellsperproc+k));
// }
// }

  // init rows we plan to use
  //for(i = k; i < k + cellsperproc; ++i)
  for(i = 0; i < cellsperproc + 2*k; ++i)
  {
    subset_v[i] = 0.0;
  }

  for(r = 1; r < rowsperproc + 1; ++r)  // +1 to offset leading ghost row
  {
    for(s = 0; s < k; ++ s)
    {
      i = (r * k) + s;

      subset_v[i] = 4 * subset_w[i];

      if(rank != 0 || r != 1)  // account for offset of ghost row
        subset_v[i] -= subset_w[i-k];
      if(s != 0)
        subset_v[i] -= subset_w[i-1];
      if(s != k-1)
        subset_v[i] -= subset_w[i+1];
      if(rank != nprocs - 1 || r != rowsperproc)  // account for offset of ghost row
        subset_v[i] -= subset_w[i+k];
    }
  }

  if(rank == 0)
  {
    MPI_Gather(subset_v + k, cellsperproc, MPI_DOUBLE, v, cellsperproc, MPI_DOUBLE,
      0, MPI_COMM_WORLD);
  }
  else
  {
    MPI_Gather(subset_v + k, cellsperproc, MPI_DOUBLE, NULL, cellsperproc,
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

