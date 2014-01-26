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

double* load_vec( char* filename, int* k );
void save_vec( int k, double* x );
double *cgsolve(double *x, int* i, double* norm, int n);
double ddot(double *v, double *w, int n);
void matvec(double *v, double *w, int n);
void daxpy(double *v, const double *w, double alpha, double beta, int n);

int rank, nprocs;

int main( int argc, char* argv[] ) {
	int writeOutX = 0;
	int n, k;
	int maxiterations = 1000;
	int niters = 0;
	double norm;
	double* b;
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

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	// Start Timer
	t1 = MPI_Wtime();
	
	// CG Solve here!
	//double x_initial[n];
	double *x=NULL;
	//x = x_initial;
if(DEBUG) printf("-calling cgsolve\n");
	x = cgsolve(x, &niters, &norm, n);
if(DEBUG) printf("-exited cgsolve\n");	
	// End Timer
	t2 = MPI_Wtime();
	
	if ( writeOutX ) {
		save_vec( k, x );
	}
		
	// Output
	printf( "Problem size (k): %d\n",k);
	if(niters > 0) {
	  printf( "Norm of the residual after %d iterations: %lf\n", niters, norm);
	}
	printf( "Elapsed time during CGSOLVE: %lf\n", t2-t1);

  if(rank == 0)
    if(x == NULL)
      printf("-x == NULL\n");
    else
	    correct = cs240_verify(x, k, t2-t1);
  else
    printf("-only rank 0 verifies\n");

	printf("Correct: %d\n", correct);
	
	// Deallocate 
	
	// if(niters > 0)
	//   free(b);

	if(rank == 0 && niters > 0)
	  free(x);
	
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
  int numrows;

if(DEBUG) printf("-started cgsolve\n");
  if(rank == 0)
  {
    x = (double *)malloc(sizeof(double) * n);
    for(i = 0; i < n; ++i)
    {
      x[i] = 0.0;
    }
  }

  numrows = n/nprocs;

if(DEBUG) printf("-initing arrays\n");
  for(i = 0; i < numrows; ++i)
  {
    r[i] = d[i] = cs240_getB(i + rank*numrows, n);
  }
 
  normb = sqrt(ddot(r, r, n));
  rtr = ddot(r, r, n);
  relres = 1.0;

if(DEBUG) printf("-starting main loop\n");
  while(relres > TARGRES && niters < MAXITERS)
  {
    ++niters;
    matvec(Ad, d, n);
if(DEBUG) printf("-left matvec\n");
    alpha = rtr / ddot(d, Ad, n);
    daxpy(x, d, 1, alpha, n);
    daxpy(r, Ad, 1, -alpha, n);
    rtrold = rtr;
    rtr = ddot(r, r, n);
    beta = rtr / rtrold;
    daxpy(d, r, beta, 1, n);
    relres = sqrt(rtr) / normb;
  }

if(DEBUG) printf("niters is %d\n", niters);

  // returning values to be printed after a run
  *iter = niters;
  *norm = relres;

  return x;
}

double ddot(double *v, double *w, int n)
{
  double prod = 0.0;
  int i;

  for(i = 0; i < n; ++i)
  {
    prod += v[i] * w[i];
  }

  return prod;
}

// Overwrites vector v with scalar1*v + scalar2*w
void daxpy(double *v, const double *w, double scalar1, double scalar2, int n)
{

  /* To parallelize - broadcast the scalars scalar1 and scalar2,
   * then split up the loop computation for the portions of v and w 
   * that you need, and combine the result at the end.
   */

  int i;
  for(i = 0; i < n; i++)
  {
    v[i] = v[i]*scalar1 + w[i]*scalar2;
  }
}

void matvec(double *v, double *w, int n)
{
  int k = (int)sqrt(n);
  int i;
  int r, s;   //row, column
  double *subset_w, *subset_v;
  int nperproc = n/nprocs;
  int kperproc = k/nprocs;
  MPI_Status status;

/*
  for(i = 0; i < n; ++i)
  {
    v[i] = 0.0;
  }
*/
  // arrays are 2k larger than we need to accomodate ghost cells from
  // neighbors
if(DEBUG) printf("-initing subsets\n");
  subset_w = (double *)malloc(sizeof(double) * (nperproc + 2*k));
  subset_v = (double *)malloc(sizeof(double) * (nperproc + 2*k));

if(DEBUG) printf("-size of subset arrays is %d\n", nperproc + 2*k);
  //double subset_v[4224];
  //double subset_w[4224];

/*
  for(i = 0; i < k; ++i)
    subset_w[i] = 0.0;
  for(i = k*(kperproc+1); i < k*(kperproc+1) + k; ++i)
    subset_w[i] = 0.0;
*/

if(DEBUG) printf("-scattering\n");
  // subset_w offset by k to give space for neighbor ghost cells
  MPI_Scatter(w, nperproc, MPI_DOUBLE, subset_w+k, nperproc, MPI_DOUBLE,
    0, MPI_COMM_WORLD);

if(DEBUG) printf("-initing ghost cells\n");
  // get ghost cells
  if(rank % 2 == 0) // even processors send then recv to up neighbor
  {
    MPI_Send(subset_w+k, k, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD);
    MPI_Recv(subset_w, k, MPI_DOUBLE, rank-1, 2, MPI_COMM_WORLD, &status);
  }
  else              // odd processors recv then send to down neighbor
  {
    MPI_Recv(subset_w+(k*(kperproc+1)), k, MPI_DOUBLE, rank+1, 1, 
      MPI_COMM_WORLD, &status);
    MPI_Send(subset_w+(k*kperproc), k, MPI_DOUBLE, rank+1, 2,
      MPI_COMM_WORLD);
  }

  if(rank % 2 == 0) // even proc send then recv to down neighbor
  {
    MPI_Send(subset_w+(k*kperproc), k, MPI_DOUBLE, rank+1, 3, 
      MPI_COMM_WORLD);
    MPI_Recv(subset_w+(k*(kperproc+1)), k, MPI_DOUBLE, rank+1, 4, 
      MPI_COMM_WORLD, &status);
  }
  else              // odd proc recv then send to up neighbor
  {
    MPI_Recv(subset_w, k, MPI_DOUBLE, rank-1, 3, MPI_COMM_WORLD, &status);
    MPI_Send(subset_w+k, k, MPI_DOUBLE, rank-1, 4, MPI_COMM_WORLD);
  }

if(DEBUG) printf("-matvec looping\n");
  for(r = 1; r < kperproc + 1; ++r)
  {
    for(s = 0; s < k; ++s)
    {
      i = (r * k) + s;
      subset_v[i] = 4 * subset_w[i];

      if(r != 0)
        subset_v[i] -= subset_w[i-k];
      if(s != 0)
        subset_v[i] -= subset_w[i-1];
      if(s != k-1)
        subset_v[i] -= subset_w[i+1];
      if(r != k-1)
        subset_v[i] -= subset_w[i+k];
    }
  }
if(DEBUG) printf("-gather\n");
  MPI_Gather(subset_v+k, nperproc, MPI_DOUBLE, v, nperproc, MPI_DOUBLE, 0,
    MPI_COMM_WORLD);

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

