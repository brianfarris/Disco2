#define TIMESTEP_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include "../Headers/Grid.h"
#include "../Headers/TimeStep.h"
#include "../Headers/header.h"

struct TimeStep *timestep_create(struct Grid * theGrid) {
  struct TimeStep *theTimeStep = (struct TimeStep *) malloc(sizeof(struct TimeStep));
  theTimeStep->t = 0.0;
  int N_r_withghost = grid_N_r(theGrid)+grid_Nghost_rmin(theGrid)+grid_Nghost_rmax(theGrid);
  int N_z_withghost = grid_N_z(theGrid)+grid_Nghost_zmin(theGrid)+grid_Nghost_zmax(theGrid);
  theTimeStep->nri = (int *)malloc(N_r_withghost*sizeof(int));
  theTimeStep->nzk = (int *)malloc(N_z_withghost*sizeof(int));
  return theTimeStep;
}

void timestep_destroy(struct TimeStep *theTimeStep) {
  free(theTimeStep->nri);
  free(theTimeStep->nzk);
  free(theTimeStep);
}


