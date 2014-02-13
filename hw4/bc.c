#include "defs.h"
#include <cilk/cilk.h>

void usage_error () {
  fprintf(stderr, "usage: ./bc <-serial|-parallel> -torus <nrows> <ncols> (computes BC for a torus)\n");
  fprintf(stderr, "usage: ./bc <-serial|-parallel> -grid <nrows> <ncols> (computes BC for a grid)\n");
  fprintf(stderr, "usage: ./bc <-serial|-parallel> (computes BC for a graph with edges from stdin\n");
  exit(-1);
}

int main(int argc, char** argv)
{
  graph* G;            // The graph data structure -- see defs.h 
  double *BC;          // BC output 
  int nrows, ncols, scale, nread, vtx, r, c, checkable;
  int *head, *tail;
  double elapsed_time;
  int torusgraph=FALSE, gridgraph=FALSE, readgraph=FALSE, serial=TRUE;

  /* --------------------*/
  /*  Parse command line */
  /* --------------------*/

  if (argc<2) usage_error();
  if (argv[1][1]=='s') serial = TRUE;
  else if (argv[1][1]=='p') serial = FALSE;
  else usage_error();
  readgraph = (argc==2);
  if (!readgraph) {
    if (argv[2][1]=='t') torusgraph = TRUE;
    else if (argv[2][1]=='g') gridgraph = TRUE;
    else usage_error();
    if (argc<4) usage_error();
    nrows = atoi(argv[3]);
    if (argc==4) ncols = nrows;
    else if (argc==5) ncols = atoi(argv[4]);
    else usage_error();
  }

  /* ------------------------------------ */
  /*  Graph construction or input         */
  /* ------------------------------------ */
  
  if (readgraph) {
    fprintf(stderr, "Reading graph edges from stdin.\n");
    nread = read_edge_list (&tail, &head);
    fprintf(stderr, "Finished reading %d graph edges.\n", nread);
    G = graph_from_edge_list (tail, head, nread);
    free(tail);
    free(head);
    fprintf(stderr, "Finished constructing graph from edge list.\n\n");
  } 
  if (gridgraph) {
    fprintf(stderr, "Generating 2D grid with %d rows and %d columns, for %d vertices in all.\n",
      nrows, ncols, nrows*ncols);
    G = (graph *) malloc(sizeof(graph));
    elapsed_time = generateGrid(G,nrows,ncols);
    fprintf(stderr, "Time to generate grid graph is %9.6lf sec.\n\n", elapsed_time);
  }
  if (torusgraph) {
    fprintf(stderr, "Generating 2D torus with %d rows and %d columns, for %d vertices in all.\n",
      nrows, ncols, nrows*ncols);
    G = (graph *) malloc(sizeof(graph));
    elapsed_time = generateTorus(G,nrows,ncols);
    fprintf(stderr, "Time to generate torus graph is %9.6lf sec.\n\n", elapsed_time);
  }
  print_CSR_graph(G);
	
  /* ------------------------------------------ */
  /*  Betweenness centrality                    */
  /* ------------------------------------------ */
  
  BC = (double *) calloc(G->nv, sizeof(double));
  if (serial) {
    fprintf(stderr, "\nRunning sequential betweenness centrality...\n");
    elapsed_time = betweennessCentrality_serial(G, BC);
  } else {
    fprintf(stderr, "\nRunning parallel betweenness centrality...\n");
    elapsed_time = betweennessCentrality_parallel(G, BC);
  }
  fprintf(stderr, "Time for betweenness centrality is %9.6f sec.\n", elapsed_time);
  fprintf(stderr, "TEPS score is %4.3e\n\n", G->nv * (G->ne) / elapsed_time );
	
  /* ------------------------------------------------------------- */
  /* Validation: Check the answer if graph is a power-of-two torus */
  /* ------------------------------------------------------------- */

  checkable = 0;
  if (torusgraph && (nrows==ncols)) {
    for (scale = 1; nrows*ncols > (1<<scale); scale++);
    checkable = nrows*ncols==(1<<scale);
  }
  if (checkable) {
			
    fprintf(stderr, "\nChecking the answer ... ");  
    double BCval = 0.5*pow(2, 3*scale/2)-pow(2, scale)+1.0;
    int failed = 0;
    for (vtx=0; vtx<G->nv; vtx++) {
      if (round(BC[vtx] - BCval) != 0) {
        failed = 1;
        break;
      }
    }
    if (failed) {
      fprintf(stderr, "sorry, answer is WRONG!\n\n");
    } else {
      fprintf(stderr, "hurray, answer is CORRECT!\n\n");
    }
  }
  else {
    fprintf(stderr, "\nAnswer checking is only available for a square, power-of-two torus.\n\n");  
  }
	
  /* ---------------------------------------------*/
  /* Print the answer if the graph is tiny enough */
  /* ---------------------------------------------*/

  if ((torusgraph||gridgraph) && ncols<=10 && nrows <= 100) {
    fprintf(stderr, "\nBetweenness centrality values:\n\n");  
    for (r = 0; r < nrows; r++) {
      for (c = 0; c < ncols; c++) {
        vtx = r*ncols+c;
        fprintf(stderr, "%6.1f ", BC[vtx]);
      }
      fprintf(stderr, "\n");
    }
    fprintf(stderr, "\n");
  } else if (G->nv < 200) {
    fprintf(stderr, "\nBetweenness centrality values:\n");  
    for (vtx = 0; vtx < G->nv; vtx++) {
      fprintf(stderr, "%6.1f ", BC[vtx]);
    }
    fprintf(stderr, "\n");
  }
	
  /* ---------*/
  /* Clean up */
  /* ---------*/

  free(BC);
  free(G->nbr);
  free(G->firstnbr);
  free(G);
  return 0;
}
