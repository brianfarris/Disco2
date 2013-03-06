#define MPISETUP_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include "../Headers/MPIsetup.h"
#include "../Headers/header.h"

struct MPIsetup *mpisetup_create(int argc,char **argv) {
  struct MPIsetup *theMPIsetup = (struct MPIsetup *) malloc(sizeof(struct MPIsetup));
  theMPIsetup->wraparound[0]=1;
  theMPIsetup->wraparound[1]=1;
  theMPIsetup->reorder = 1;
  MPI_Init(&argc,&argv);
  return theMPIsetup;
}

void mpisetup_destroy(struct MPIsetup *theMPIsetup) {
  MPI_Barrier(grid_comm);
  MPI_Finalize();
  free(theMPIsetup);
}


