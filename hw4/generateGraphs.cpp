#include "defs.h"

double generateGrid(graph* G, int nrows, int ncols) {
  //  Generate a grid with nrows*ncols vertices
  //  Returns elapsed time in seconds
  double elapsed_time;
  int r, c, k, vtx, nv, ne;

  elapsed_time = get_seconds();
  nv = nrows*ncols;
  ne = 4*nv - 2*(nrows+ncols);  
  G->nv = nv;
  G->ne = ne;                    
  G->nbr = (int *) malloc(ne*sizeof(int));
  G->firstnbr = (int *) malloc((nv+1)*sizeof(int));
  k = 0;
  for (r = 0; r < nrows; r++) {
    for (c = 0; c < ncols; c++) {
      vtx = r*ncols+c;
      G->firstnbr[vtx] = k;
      if (r!=0)       G->nbr[k++] = vtx-ncols; // all but first grid row
      if (c!=0)       G->nbr[k++] = vtx-1;     // all but first grid col
      if (c!=ncols-1) G->nbr[k++] = vtx+1;     // all but last grid col
      if (r!=nrows-1) G->nbr[k++] = vtx+ncols; // all but last grid row
    }
  }
  assert(k==ne);
  G->firstnbr[nv] = ne;
  elapsed_time = get_seconds() - elapsed_time; 
  return elapsed_time;
}


double generateTorus(graph* G, int nrows, int ncols) {
  //  Generate a torus with nrows*ncols vertices
  //  Returns elapsed time in seconds
  double elapsed_time;
  int r, c, k, vtx, nv, ne;

  elapsed_time = get_seconds();
  nv = nrows*ncols;
  ne = 4*nv;  
  G->nv = nv;
  G->ne = ne;                    
  G->nbr = (int *) malloc(ne*sizeof(int));
  G->firstnbr = (int *) malloc((nv+1)*sizeof(int));
  k = 0;
  for (r = 0; r < nrows; r++) {
    for (c = 0; c < ncols; c++) {
      vtx = r*ncols+c;
      G->firstnbr[vtx] = k;
      if (r==0)       G->nbr[k++] = vtx-ncols+nv;  // first grid row
        else          G->nbr[k++] = vtx-ncols; 
      if (c==0)       G->nbr[k++] = vtx-1+ncols;   // first grid col
        else          G->nbr[k++] = vtx-1;     
      if (c==ncols-1) G->nbr[k++] = vtx+1-ncols;   // last grid col
        else          G->nbr[k++] = vtx+1;
      if (r==nrows-1) G->nbr[k++] = vtx+ncols-nv;  // last grid row
        else          G->nbr[k++] = vtx+ncols; 
    }
  }
  assert(k==ne);
  G->firstnbr[nv] = ne;
  elapsed_time = get_seconds() - elapsed_time; 
  return elapsed_time;
}
