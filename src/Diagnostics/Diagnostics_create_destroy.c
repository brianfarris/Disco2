#define DIAGNOSTICS_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include "../Headers/Grid.h"
#include "../Headers/Diagnostics.h"
#include "../Headers/TimeStep.h"
#include "../Headers/header.h"

struct Diagnostics *diagnostics_create(struct Grid * theGrid, struct TimeStep * theTimeStep) {
  struct Diagnostics * theDiagnostics = (struct Diagnostics *) malloc(sizeof(struct Diagnostics));
  
  theDiagnostics->VectorDiag = malloc(sizeof(double *) * grid_N_r_global(theGrid));
  theDiagnostics->VectorDiag[0] = malloc(sizeof(double) * grid_N_r_global(theGrid) * theDiagnostics->NUM_DIAG);
  int i,n;
  for (i=0;i<grid_N_r_global(theGrid);++i){
    theDiagnostics->VectorDiag[i] = theDiagnostics->VectorDiag[0]+i*theDiagnostics->NUM_DIAG;
  }

  theDiagnostics->ScalarDiag = malloc(sizeof(double)*theDiagnostics->NUM_DIAG);

  for (i=0;i<grid_N_r_global(theGrid);++i){
    for (n=0;n<theDiagnostics->NUM_DIAG;++n){
       theDiagnostics->VectorDiag[i][n]=0.0;
    }
  }
  for (n=0;n<theDiagnostics->NUM_DIAG;++n){
       theDiagnostics->ScalarDiag[n]=0.0;
    }

  theDiagnostics->toutprev = timestep_get_t(theTimeStep); 
  theDiagnostics->toutprev_dump = timestep_get_t(theTimeStep); 
  return theDiagnostics;
}

void diagnostics_destroy(struct Diagnostics *theDiagnostics,struct Grid * theGrid) {
  int num_r_points_global = grid_N_r_global(theGrid);
  int i,n;
  free(theDiagnostics->ScalarDiag);
  free(theDiagnostics->VectorDiag[0]);
  free(theDiagnostics->VectorDiag);
  free(theDiagnostics);
}


