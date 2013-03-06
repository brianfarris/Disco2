#define GRID_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include "../Headers/Grid.h"
#include "../Headers/MPIsetup.h"
#include "../Headers/header.h"


struct Grid *grid_create(struct MPIsetup * theMPIsetup) {
  struct Grid *theGrid = (struct Grid *) malloc(sizeof(struct Grid));
  return theGrid;
}

void grid_destroy(struct Grid *theGrid) {
  free(theGrid->r_faces);
  free(theGrid->z_faces);
  free(theGrid->N_p);
  free(theGrid);
}


