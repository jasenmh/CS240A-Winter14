#include "defs.h"
#include <cilk/cilk.h>
#include <cilk/reducer_opadd.h>
#include <strings.h>

#define MAX_THREADS 320
#define DEBUG 1

int numvertsingraph;  // need to make this a global variable so reducer identity function can handle varying sized graphs

void centrality_update(double *left, int node, double cu)
{
  left[node] += cu;
}

void identity_wrapper(void *reducer, void *bcr)
{
  double *tmp = (double *)malloc(sizeof(double) * numvertsingraph);
  bzero(tmp, numvertsingraph * sizeof(double));
  bcr = (void *)tmp;
}

void reducer_wrapper(void *reducer, void *left, void *right)
{
  int i;
  double *dl = (double *)left;
  double *dr = (double *)right;

  for(i = 0; i < numvertsingraph; ++i)
  {
    dl[i] += dr[i];
  }

}

void destroy_wrapper(void *reducer, void *p)
{
  if(p == NULL)
    return;

  //free(p);
}

double betweennessCentrality_parallel(graph* G, double* BC) {
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
  int continueforever;

  // create and initialize our custom reducer
if(DEBUG) printf("- creating and initializing reducer\n");
  numvertsingraph = G->nv;
  CILK_C_DECLARE_REDUCER(double *) my_bcr =
    CILK_C_INIT_REDUCER(double *,
      reducer_wrapper,
      identity_wrapper,
//      destroy_wrapper);
      __cilkrts_hyperobject_noop_destroy
      );
  //BC_initialize(&REDUCER_VIEW(my_bcr));
  bzero(&REDUCER_VIEW(my_bcr), numvertsingraph * sizeof(double));

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
  if(in_degree != NULL)
    free(in_degree);
  if(numEdges != NULL)
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

  // register our custom reducer
if(DEBUG) printf("- registering reducer\n");
  CILK_C_REGISTER_REDUCER(my_bcr);

  continueforever = 0;

  // cilk_for this for loop
  //for (p=0; p<n; p++) {
  cilk_for(p = 0; p < n; ++p)
  {

		i = Srcs[p];
		if (G->firstnbr[i+1] - G->firstnbr[i] == 0) {
			continue;
		} else {
			num_traversals++;
		}

		if ((num_traversals == numV + 1) || continueforever) {
      // can't break out of a cilk_for loop
			//break;
			continueforever = 1;
      continue;
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
        // replace this BC with our reducer
				//BC[w] += del[w];
//if(DEBUG) printf("- updating centrality in reducer\n");
				//BC_centrality_update(&REDUCER_VIEW(my_bcr), w, del[w]);
				//*(REDUCER_VIEW(my_bcr) + w) = *(REDUCER_VIEW(my_bcr) + w) + del[w];
        centrality_update(REDUCER_VIEW(my_bcr), w, del[w]);
			}

			phase_num--;
		}
		
		for (j=0; j<count; j++) {
			w = S[j];
			d[w] = -1;
			del[w] = 0;
			P[w].count = 0;
		}
if(DEBUG) printf("- about to exit loop\n");
  }
if(DEBUG) printf("- exited loop\n");
  /***********************************/
  /*** END OF MAIN LOOP **************/
  /***********************************/

  // unregister reducer
if(DEBUG) printf("- unregistering reducer\n");
  CILK_C_UNREGISTER_REDUCER(my_bcr);

  // copy the values from our reducer into BC
if(DEBUG) printf("- exporting betweennessess\n");
	//BC_export_betweennessess(&REDUCER_VIEW(my_bcr), BC);
	for(i = 0; i < numvertsingraph; ++i)
  {
    BC[i] = REDUCER_VIEW(my_bcr)[i];
  }

  if(S != NULL)
    free(S);
  if(pListMem != NULL)
    free(pListMem);
  if(P != NULL)
    free(P);
  if(sig != NULL)
    free(sig);
  if(d != NULL)
    free(d);
  if(del != NULL)
    free(del);
  if(start != NULL)
    free(start);
  if(end != NULL)
    free(end);

  elapsed_time = get_seconds() - elapsed_time;

  if(Srcs != NULL)
    free(Srcs);

  // TODO: profit!
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
