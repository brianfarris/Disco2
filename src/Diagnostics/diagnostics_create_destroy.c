#define DIAGNOSTICS_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include "../Headers/Sim.h"
#include "../Headers/Diagnostics.h"
#include "../Headers/TimeStep.h"
#include "../Headers/header.h"

struct Diagnostics *diagnostics_create(struct Sim * theSim, struct TimeStep * theTimeStep) {
  struct Diagnostics * theDiagnostics = (struct Diagnostics *) malloc(sizeof(struct Diagnostics));

  theDiagnostics->NUM_DIAG = 2;

  theDiagnostics->VectorDiag = malloc(sizeof(double *) * sim_N_r_global(theSim));
  theDiagnostics->VectorDiag[0] = malloc(sizeof(double) * sim_N_r_global(theSim) * (theDiagnostics->NUM_DIAG+1));
  int i,n;
  for (i=0;i<sim_N_r_global(theSim);++i){
    theDiagnostics->VectorDiag[i] = theDiagnostics->VectorDiag[0]+i*(theDiagnostics->NUM_DIAG+1);
  }

  theDiagnostics->ScalarDiag = malloc(sizeof(double)*theDiagnostics->NUM_DIAG);

  for (i=0;i<sim_N_r_global(theSim);++i){
    for (n=0;n<(1+theDiagnostics->NUM_DIAG);++n){
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

void diagnostics_destroy(struct Diagnostics *theDiagnostics,struct Sim * theSim) {
  int num_r_points_global = sim_N_r_global(theSim);
  int i,n;
  free(theDiagnostics->ScalarDiag);
  free(theDiagnostics->VectorDiag[0]);
  free(theDiagnostics->VectorDiag);
  free(theDiagnostics);
}


