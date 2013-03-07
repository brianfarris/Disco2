#define TIMESTEP_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include "../Headers/Grid.h"
#include "../Headers/TimeStep.h"
#include "../Headers/header.h"

struct TimeStep *timestep_create(struct Grid * theGrid) {
  struct TimeStep *theTimeStep = (struct TimeStep *) malloc(sizeof(struct TimeStep));
  theTimeStep->t = 0.0;
  theTimeStep->nri = (int *)malloc(grid_N_r(theGrid)*sizeof(int));
  theTimeStep->nzk = (int *)malloc(grid_N_z(theGrid)*sizeof(int));
  return theTimeStep;
}

void timestep_destroy(struct TimeStep *theTimeStep) {
  free(theTimeStep->nri);
  free(theTimeStep->nzk);
  free(theTimeStep);
}


