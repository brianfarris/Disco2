#define SIM_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include "../Headers/Sim.h"
#include "../Headers/MPIsetup.h"
#include "../Headers/header.h"


struct Sim *sim_create(struct MPIsetup * theMPIsetup) {
  struct Sim *theSim = (struct Sim *) malloc(sizeof(struct Sim));
  return theSim;
}

void sim_destroy(struct Sim *theSim) {
  free(theSim->r_faces);
  free(theSim->z_faces);
  free(theSim->N_p);
  free(theSim);
}


