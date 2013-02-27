#define TIMESTEP_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include "../Headers/TimeStep.h"
#include "../Headers/header.h"

struct TimeStep *timestep_create() {
  struct TimeStep *theTimeStep = (struct TimeStep *) malloc(sizeof(struct TimeStep));
  theTimeStep->t = 0.0;
  theTimeStep->T_MAX = 50.0;
  theTimeStep->NUM_CHECKPOINTS = 100;
  return theTimeStep;
}

void timestep_destroy(struct TimeStep *theTimeStep) {
  free(theTimeStep);
}


