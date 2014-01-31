#include <cilk/cilk.h>
#include <cilk/reducer_opadd.h>

#include <algorithm>
#include <numeric>
#include <string>
#include <iostream>
#include <cmath>
#include <iterator>
#include <functional>

#include "example_util_gettime.h"

#define COARSENESS 100
#define ITERS 10

double rec_cilkified(double * a, double * b, int n)
{
  int split_on;
  int i;
  double suma, sumb;

  if(n <= COURSENESS) // go serial
  {
    suma = 0.0;

    for(i = 0; i < n; ++i)
    {
      suma += a[i] * b[i];
    }

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
  int outer, inner;
  int npercourseness = n / COURSENESS;
  int inneridx;
  double sum = 0.0;
  double prods[n];

  // schedule segments of the vectors to be multiplied in parallel
  cilk_for(outer = 0; outer < npercourseness; ++outer)
  {
    for(inner = 0; inner < COURSENESS; ++inner)
    {
      inneridx = (outer * COURSENESS) + inner;

      prods[inneridx] = a[inneridx] * b[inneridx];
    }
  }

  // sum the products in serial
  for(outer = 0; outer < n; ++outer)
  {
    sum += prods[outer];
  }

	return sum;
}


double hyperobject_cilkified(double * a, double * b, int n)
{
  cilk::hyperobject<cilk::reducer_opadd<double> > sum;
  int i;

  cilk_for(i = 0; i < n; ++i)
  {
    sum += a[i] * b[i];
  }

  return sum.get_value();
}


int close(double x, double y, int n)
{
        double relative_diff = (x>y? (x-y)/x : (y-x)/x);
        return (relative_diff < sqrt((double) n) * exp2(-42))? 1 : 0;
}


// A simple test harness 
int inn_prod_driver(int n)
{
	double * a = new double[n];
	double * b = new double[n];
	for (int i = 0; i < n; ++i)
	{
        	a[i] = i;
		b[i] = i;
	}
    	std::random_shuffle(a, a + n);
	std::random_shuffle(b, b + n);

	double seqresult = std::inner_product(a, a+n, b, 0);	

	long t1 = example_get_time();
	for(int i=0; i< ITERS; ++i)
	{
		seqresult = std::inner_product(a, a+n, b, 0);	
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
