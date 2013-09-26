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

void timestep_set_dt(struct TimeStep * theTimeStep, struct Cell *** theCells, struct Sim * theSim, struct GravMass * theGravMasses){
  //theTimeStep->dt = (1.95/(1+exp(-theTimeStep->t/0.01))-0.95) * cell_mindt(theCells,theSim);
  theTimeStep->dt = cell_mindt(theCells,theSim, theGravMasses);
  if( theTimeStep->t+theTimeStep->dt > sim_get_T_MAX(theSim) ) {
    theTimeStep->dt = sim_get_T_MAX(theSim)-theTimeStep->t;
  }
}
void timestep_update_t(struct TimeStep * theTimeStep){
  theTimeStep->t += theTimeStep->dt;
  time_global = theTimeStep->t;
}
void timestep_set_RK(struct TimeStep * theTimeStep,double RK){
  theTimeStep->RK = RK;
}
void timestep_set_t(struct TimeStep * theTimeStep,double t){
  theTimeStep->t = t;
}
void timestep_set_nri(struct TimeStep * theTimeStep,int i,int n){
  theTimeStep->nri[i] = n;
}
void timestep_set_nzk(struct TimeStep * theTimeStep,int k,int n){
  theTimeStep->nzk[k] = n;
}
