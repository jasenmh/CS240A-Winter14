#include "defs.h"
#include <cilk/cilk.h>
#include <cilk/reducer_opadd.h>

#define MAX_THREADS 320

double betweennessCentrality_parallel(graph* G, double* BC) {
  double elapsed_time;
  int i;
  double c;
  cilk_for(i=0; i<10;i++){
      c+= i;
  }
  printf("c = %f\n",c);
  return elapsed_time;
}

/*
 * Serial Version
 *
 */
double betweennessCentrality_serial(graph* G, double* BC) {
  int *S; 	/* stack of vertices in order of distance from s. Also, implicitly, the BFS queue */
  plist* P;  	/* predecessors of vertex v on shortest paths from s */
  double* sig; 	/* No. of shortest paths */
  int* d; 	/* Length of the shortest path between every pair */
  double* del; 	/* dependency of vertices */
  int *in_degree, *numEdges;
  int *pListMem;	
  int* Srcs; 
  int *start, *end;
  int seed = 2387;
  double elapsed_time;
  int i, j, k, p, count, myCount;
  int v, w, vert;
  int numV, num_traversals, n, m, phase_num;

  /* numV: no. of vertices to run BFS from = 2^K4approx */
  //numV = 1<<K4approx;
  n = G->nv;
  m = G->ne;
  numV = n;

  /* Permute vertices */
  Srcs = (int *) malloc(n*sizeof(int));
  for (i=0; i<n; i++) {
    Srcs[i] = i;
  }

  /* Start timing code from here */
  elapsed_time = get_seconds();

  /* Initialize predecessor lists */
  /* Number of predecessors of a vertex is at most its in-degree. */
  P = (plist *) calloc(n, sizeof(plist));
  in_degree = (int *) calloc(n+1, sizeof(int));
  numEdges = (int *) malloc((n+1)*sizeof(int));
  for (i=0; i<m; i++) {
    v = G->nbr[i];
    in_degree[v]++;
  }
  prefix_sums(in_degree, numEdges, n);
  pListMem = (int *) malloc(m*sizeof(int));
  for (i=0; i<n; i++) {
    P[i].list = pListMem + numEdges[i];
    P[i].degree = in_degree[i];
    P[i].count = 0;
  }
  free(in_degree);
  free(numEdges);
	
  /* Allocate shared memory */ 
  S   = (int *) malloc(n*sizeof(int));
  sig = (double *) malloc(n*sizeof(double));
  d   = (int *) malloc(n*sizeof(int));
  del = (double *) calloc(n, sizeof(double));
	
  start = (int *) malloc(n*sizeof(int));
  end = (int *) malloc(n*sizeof(int));

  num_traversals = 0;
  myCount = 0;

  for (i=0; i<n; i++) {
    d[i] = -1;
  }
	
  /***********************************/
  /*** MAIN LOOP *********************/
  /***********************************/
  for (p=0; p<n; p++) {

		i = Srcs[p];
		if (G->firstnbr[i+1] - G->firstnbr[i] == 0) {
			continue;
		} else {
			num_traversals++;
		}

		if (num_traversals == numV + 1) {
			break;
		}
		
		sig[i] = 1;
		d[i] = 0;
		S[0] = i;
		start[0] = 0;
		end[0] = 1;
		
		count = 1;
		phase_num = 0;

		while (end[phase_num] - start[phase_num] > 0) {
				myCount = 0;
				// BFS to destination, calculate distances, 
				int vert;
				for ( vert = start[phase_num]; vert < end[phase_num]; vert++ ) {
					v = S[vert];
					int j;
					for ( j=G->firstnbr[v]; j<G->firstnbr[v+1]; j++ ) {
						w = G->nbr[j];
						if (v != w) {

							/* w found for the first time? */ 
							if (d[w] == -1) {
								//printf("n=%d, j=%d, start=%d, end=%d, count=%d, vert=%d, w=%d, v=%d\n",n,j,start[phase_num],end[phase_num],myCount,vert,w,v);
								S[end[phase_num] + myCount] = w;
								myCount++;
								d[w] = d[v] + 1; 
								sig[w] = sig[v]; 
								P[w].list[P[w].count++] = v;
							} else if (d[w] == d[v] + 1) {
								sig[w] += sig[v]; 
								P[w].list[P[w].count++] = v;
							}
						
						}
					}
	 			}
			
				/* Merge all local stacks for next iteration */
				phase_num++; 
				
				start[phase_num] = end[phase_num-1];
				end[phase_num] = start[phase_num] + myCount;
			
				count = end[phase_num];
		}
 	
		phase_num--;

		while (phase_num > 0) {
			for (j=start[phase_num]; j<end[phase_num]; j++) {
				w = S[j];
				for (k = 0; k < P[w].count; k++) {
					v = P[w].list[k];
					del[v] = del[v] + sig[v]*(1+del[w])/sig[w];
				}
				BC[w] += del[w];
			}

			phase_num--;
		}
		
		for (j=0; j<count; j++) {
			w = S[j];
			d[w] = -1;
			del[w] = 0;
			P[w].count = 0;
		}
  }
  /***********************************/
  /*** END OF MAIN LOOP **************/
  /***********************************/
 

	
  free(S);
  free(pListMem);
  free(P);
  free(sig);
  free(d);
  free(del);
  free(start);
  free(end);
  elapsed_time = get_seconds() - elapsed_time;
  free(Srcs);

  return elapsed_time;
}
