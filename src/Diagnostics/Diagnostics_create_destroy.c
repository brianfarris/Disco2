#define DIAGNOSTICS_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include "../Headers/Grid.h"
#include "../Headers/Diagnostics.h"
#include "../Headers/TimeStep.h"
#include "../Headers/header.h"

struct Diagnostics *diagnostics_create(struct Grid * theGrid, struct TimeStep * theTimeStep) {
  int NUM_DIAG = 2;
  struct Diagnostics * theDiagnostics = (struct Diagnostics *) malloc(sizeof(struct Diagnostics));
  
  theDiagnostics->VectorDiag = malloc(sizeof(double *) * grid_N_r_global(theGrid));
  theDiagnostics->VectorDiag[0] = malloc(sizeof(double) * grid_N_r_global(theGrid) * NUM_DIAG);
  int i;
  for (i=0;i<grid_N_r_global(theGrid);++i){
    theDiagnostics->VectorDiag[i] = theDiagnostics->VectorDiag[0]+i*NUM_DIAG;
  }

  theDiagnostics->ScalarDiag = malloc(sizeof(double)*NUM_DIAG);

  theDiagnostics->toutprev = timestep_get_t(theTimeStep); 
  return theDiagnostics;
}

void diagnostics_destroy(struct Diagnostics *theDiagnostics,struct Grid * theGrid) {
  int NUM_DIAG = 2;
  int num_r_points_global = grid_N_r_global(theGrid);
  //printf("theDiagnostics->ScalarDiag[0]: %e\n",theDiagnostics->ScalarDiag[0]);
  int i,n;
  /*
  for (i=0;i<num_r_points_global;++i){
    printf("theDiagnostics->VectorDiag[%d][NUM_DIAG]: %e\n",i,theDiagnostics->VectorDiag[i][0]);
  }
  */
  free(theDiagnostics->ScalarDiag);
  free(theDiagnostics->VectorDiag[0]);
  free(theDiagnostics->VectorDiag);
  free(theDiagnostics);
}


