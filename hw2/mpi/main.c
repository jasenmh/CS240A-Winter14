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

double* load_vec( char* filename, int* k );
void save_vec( int k, double* x );
double *cgsolve(int n);
double ddot(double *v, double *w, int n);
double *daxpy();
void matvec(double *v, const double *w, int n);

int main( int argc, char* argv[] ) {
	int writeOutX = 0;
	int n, k;
	int iterations = 1000;
	double norm;
	double* b;
	double* x;
	double time;
	double t1, t2;
	
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
	
	// End Timer
	t2 = MPI_Wtime();
	
	if ( writeOutX ) {
		save_vec( k, x );
	}
		
	// Output
	printf( "Problem size (k): %d\n",k);
	printf( "Norm of the residual after %d iterations: %lf\n",iterations,norm);
	printf( "Elapsed time during CGSOLVE: %lf\n", t1-t2);
	
	// Deallocate 
	free(b);
	free(x);
	
	MPI_Finalize();
	
	return 0;
}

/*
 * Team functions
 *
 */

double *cgsolve(int n)
{
  int niters = 0;   // number of iterations
  int possmax = 5*sqrt(n);  // possible max. iterations
  int MAXITERS = (1000 > possmax) ? 1000 : possmax;
  double TARGRES = 1.0e-6;  // target residual
  double relres;  // relative residual
  double x[n];
  double b[n];
  double r[n];
  double d[n];  // direction
  double Ad[n]
  double alpha, beta;
  double rtr;
  double rtrold;
  double normb;
  int i;

  for(i = 0; i < n; ++i)
  {
    x[i] = 0;
    b[i] = r[i] = d[i] = getb(i);
  }
 
  rtr = ddot(r, r);
  normb = sqrt(ddot(b, b));
  relres = 1;

  while(relres < TARGRES && niter < MAXITERS)
  {
    ++niters;
    matvec(Ad, d, n);
    alpha = rtr / ddot(d, Ad);
    x = daxpy(x, d, alpha, n);
    r = daxpy(r, Ad, -alpha, n);
    rtrold = rtr;
    rtr = ddot(r, r);
    beta = rtr / rtrold;
    relres = sqrt(rtr) / normb;
  }

}

double ddot(double *v, double *w, int n)
{
  double prod = 0;
  int i;

  for(i = 0; i < n; ++i)
  {
    p += v[i] * w[i];
  }

  return p;
}

void matvec(double *v, const double *w, int n)
{
  int k = sqrt(n);
  int i;
  int r, s;   //row, column

/* id don't think i need to init this
  for(i = 0; i < n; ++i)
  {
    v[i] = 0.0;
  }
*/

  for(r = 0; r < k; ++r)
  {
    for(s = 0; s < k; ++s)
    {
      i = ((r - 1) * k) + s;
      v[i] = 4 * w[i];

      if(r != 1)
        v[i] -= w[i-k];
      if(s != 1)
        v[i] -= w[i-1];
      if(s != k)
        v[i] -= w[i+1];
      if(r != k)
        v[i] -= w[i+k];
    }
  }

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

