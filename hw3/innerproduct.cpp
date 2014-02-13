// Jasen Hall and Kyle Jorgensen, Assignment 3 for CS240A, 1 Feb 2014

#include <cilk/cilk.h>
#include <cilk/reducer_opadd.h>
#include <cilk/cilk_api.h>

#include <algorithm>
#include <numeric>
#include <string>
#include <iostream>
#include <cmath>
#include <iterator>
#include <functional>

#include "example_util_gettime.h"

#define COARSENESS 1000
#define ITERS 10

double rec_cilkified(double * a, double * b, int n)
{
  int split_on;
  int i;
  double suma, sumb;

  if(n <= COARSENESS) // go serial
  {
    suma = std::inner_product(a, a+n, b, (double)0);
    return suma;
  }
  else  // split and recurse
  {
    split_on = n / 2; // truncting decimal
    suma = cilk_spawn rec_cilkified(a, b, split_on);
    sumb = rec_cilkified(a + split_on, b + split_on, n - split_on);
    cilk_sync;

    return suma + sumb;
  }

	return 0;
}

double loop_cilkified(double * a, double * b, int n)
{
  int npercoarseness = n / COARSENESS;
  int remainder = n % COARSENESS;
  int inneridx;
  double sum = 0.0;
  // create a new array initialized with zeros
  double * partialProd = new double[npercoarseness]();
  double * extraPartialProd;

  if (remainder != 0) // if n is not a factor of the coarseness, we need an extra array
  	extraPartialProd = new double[remainder]();

  // schedule segments of the vectors to be multiplied in parallel
  cilk_for(int outer = 0; outer < npercoarseness; ++outer)
  {
    for(int inner = 0; inner < COARSENESS; ++inner)
    {
      inneridx = (outer * COARSENESS) + inner;

      partialProd[outer] += a[inneridx] * b[inneridx];
    }
  }

  if(remainder != 0) // compute the remainder of the arrays
    cilk_for(int i = 1; i <= remainder; i++)
      { extraPartialProd[i-1] += a[n-i] * b[n-i]; }

  // sum the products in serial
  for(int outer = 0; outer < npercoarseness; ++outer)
  {
    sum += partialProd[outer];
  }

  for(int i = 0; i < remainder; ++i)
  {
  	sum += extraPartialProd[i];
  }

	return sum;
}


double hyperobject_cilkified(double * a, double * b, int n)
{
  cilk::reducer_opadd<double> parallel_sum;

  cilk_for(int i = 0; i < n; ++i)
  {
    parallel_sum += a[i] * b[i];
  }

  return parallel_sum.get_value();
}


int close(double x, double y, int n)
{
        double relative_diff = (x>y? (x-y)/x : (y-x)/x);
        return (relative_diff < sqrt((double) n) * exp2(-42))? 1 : 0;
}


// A simple test harness 
int inn_prod_driver(int n)
{
  __cilkrts_set_param("nworkers","4");

	double * a = new double[n];
	double * b = new double[n];
	for (int i = 0; i < n; ++i)
	{
        	a[i] = i;
		b[i] = i;
	}
    	std::random_shuffle(a, a + n);
	std::random_shuffle(b, b + n);

	double seqresult = std::inner_product(a, a+n, b, (double)0);	

	long t1 = example_get_time();
	for(int i=0; i< ITERS; ++i)
	{
		seqresult = std::inner_product(a, a+n, b, (double)0);	
	}
	long t2 = example_get_time();

	double seqtime = (t2-t1)/(ITERS*1000.f);
	std::cout << "Sequential time: " << seqtime << " seconds" << std::endl;	
	
	/***********************************************************/
	/********  START TESTING RECURSIVE CILKFIED VERSION  *******/
	/***********************************************************/

	double parresult = rec_cilkified(a, b, n);   
	t1 = example_get_time();
	for(int i=0; i< ITERS; ++i)
	{
		parresult = rec_cilkified(a, b, n);   
	}
 	t2 = example_get_time();

	double partime = (t2-t1)/(ITERS*1000.f);
	std::cout << "Recursive cilkified time:" << partime << " seconds" << std::endl;
	std::cout << "Speedup is: " << seqtime/partime << std::endl;
	std::cout << "Sequential result is: "<<seqresult<<std::endl;
	std::cout << "Recursive cilkified result is: "<<parresult<<std::endl;
	std::cout << "Result is " << (close(seqresult,parresult,n)  ? "correct":"incorrect") << std::endl; 
	
	/****************************************************************/
	/********  START TESTING NESTED LOOPED CILKIFIED VERSION  *******/
	/****************************************************************/
	parresult = loop_cilkified(a, b, n);   
	
	t1 = example_get_time();
	for(int i=0; i< ITERS; ++i)
	{
		//parresult = loop_cilkified(a, b, n);   
 	        parresult = loop_cilkified(a, b, n);   
	}
 	t2 = example_get_time();


	partime = (t2-t1)/(ITERS*1000.f);
	std::cout << "Nested loop cilkified time: " << partime << " seconds" << std::endl;
	std::cout << "Speedup is: " << seqtime/partime << std::endl;
	std::cout << "Sequential result is: "<<seqresult<<std::endl;
	std::cout << "Loop cilkified result is: "<<parresult<<std::endl;
	std::cout << "Result is " << (close(seqresult,parresult,n)  ? "correct":"incorrect") << std::endl; 
	
	/**************************************************************/
	/********  START TESTING HYPEROBJECT CILKIFIED VERSION  *******/
	/**************************************************************/

	parresult = hyperobject_cilkified(a, b, n);   
	
	t1 = example_get_time();
	for(int i=0; i< ITERS; ++i)
	{
		parresult = hyperobject_cilkified(a, b, n);   
	}
 	t2 = example_get_time();

	partime = (t2-t1)/(ITERS*1000.f);
	std::cout << "Hyperobject cilkified time:" << partime << " seconds" << std::endl;
	std::cout << "Speedup is: " << seqtime/partime << std::endl;
	std::cout << "Sequential result is: "<<seqresult<<std::endl;
	std::cout << "Hyperobject result is: "<<parresult<<std::endl;
	std::cout << "Result is " << (close(seqresult,parresult,n)  ? "correct":"incorrect") << std::endl; 
    	
        
	delete [] a;
	delete [] b;
    	return 0;
}

int main(int argc, char* argv[])
{
    int n = 1 * 1000 * 1000;
    if (argc > 1) {
        n = std::atoi(argv[1]);
    }

    return inn_prod_driver(n);
}
