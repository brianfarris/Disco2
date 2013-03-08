#define DIAGNOSTICS_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include "../Headers/Grid.h"
#include "../Headers/Diagnostics.h"
#include "../Headers/header.h"

struct Diagnostics *diagnostics_create(struct Grid * theGrid) {
  int NUM_DIAG = 2;
  int N_r_noghost = grid_N_r(theGrid)-grid_Nghost_rmin(theGrid)-grid_Nghost_rmax(theGrid);
  struct Diagnostics * theDiagnostics = (struct Diagnostics *) malloc(sizeof(struct Diagnostics));
  
  theDiagnostics->VectorDiag = malloc(sizeof(double *) * N_r_noghost);
  theDiagnostics->VectorDiag[0] = malloc(sizeof(double) * N_r_noghost * NUM_DIAG);
  int i;
  for (i=0;i<N_r_noghost;++i){
    theDiagnostics->VectorDiag[i] = theDiagnostics->VectorDiag[0]+i*NUM_DIAG;
  }

  theDiagnostics->ScalarDiag = malloc(sizeof(double)*NUM_DIAG);

  return theDiagnostics;
}

void diagnostics_destroy(struct Diagnostics *theDiagnostics) {
  free(theDiagnostics->ScalarDiag);
  free(theDiagnostics->VectorDiag[0]);
  free(theDiagnostics->VectorDiag);
  free(theDiagnostics);
}


