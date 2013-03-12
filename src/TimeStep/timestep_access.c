#define TIMESTEP_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/TimeStep.h"
#include "../Headers/Riemann.h"
#include "../Headers/Sim.h"
#include "../Headers/Cell.h"
#include "../Headers/GravMass.h"
#include "../Headers/Face.h"
#include "../Headers/TimeStep.h"
#include "../Headers/MPIsetup.h"
#include "../Headers/header.h"

double timestep_get_t(struct TimeStep * theTimeStep){
  return(theTimeStep->t);
}
double timestep_dt(struct TimeStep * theTimeStep){
  return(theTimeStep->dt);
}
/*
int timestep_nri(struct TimeStep * theTimeStep, int i){
  return(theTimeStep->nri[i]);
}
int timestep_nzk(struct TimeStep * theTimeStep, int k){
  return(theTimeStep->nzk[k]);
}
*/
int timestep_n(struct TimeStep * theTimeStep, int index, int direction){
  if (direction==R_DIR){
    return(theTimeStep->nri[index]);
  } else if (direction==Z_DIR){
    return(theTimeStep->nzk[index]);
  } else {
    printf("ERROR\n");
    exit(0);
  }
}
int timestep_Nfr(struct TimeStep * theTimeStep){
  return(theTimeStep->Nfr);
}
int timestep_Nfz(struct TimeStep * theTimeStep){
  return(theTimeStep->Nfz);
}


