#define TIMESTEP_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include "../Headers/Sim.h"
#include "../Headers/TimeStep.h"
#include "../Headers/header.h"

struct TimeStep *timestep_create(struct Sim * theSim) {
  struct TimeStep *theTimeStep = (struct TimeStep *) malloc(sizeof(struct TimeStep));
  theTimeStep->t = 0.0;
  theTimeStep->nri = (int *)malloc(sim_N(theSim,R_DIR)*sizeof(int));
  theTimeStep->nzk = (int *)malloc(sim_N(theSim,Z_DIR)*sizeof(int));
  return theTimeStep;
}

void timestep_destroy(struct TimeStep *theTimeStep) {
  free(theTimeStep->nri);
  free(theTimeStep->nzk);
  free(theTimeStep);
}


