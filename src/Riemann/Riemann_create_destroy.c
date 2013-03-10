#define RIEMANN_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Sim.h"
#include "../Headers/Riemann.h"
#include "../Headers/header.h"

struct Riemann *riemann_create(struct Sim *theSim){
  struct Riemann * theRiemann = malloc(sizeof(struct Riemann));
  theRiemann->primL = malloc(sizeof(double)*sim_NUM_Q(theSim));
  theRiemann->primR = malloc(sizeof(double)*sim_NUM_Q(theSim));
  theRiemann->Uk = malloc(sizeof(double)*sim_NUM_Q(theSim));
  theRiemann->Ustar = malloc(sizeof(double)*sim_NUM_Q(theSim));
  theRiemann->F = malloc(sizeof(double)*sim_NUM_Q(theSim));
  int q;
  for (q=0;q<sim_NUM_Q(theSim);++q){
    theRiemann->primL[q]=0.;
    theRiemann->primR[q]=0.;
    theRiemann->Uk[q]=0.;
    theRiemann->Ustar[q]=0.;
    theRiemann->F[q]=0.;
  }
  return(theRiemann);
}

void riemann_destroy(struct Riemann * theRiemann){
  free(theRiemann->F);
  free(theRiemann->primL);
  free(theRiemann->primR);
  free(theRiemann->Uk);
  free(theRiemann->Ustar);
  free(theRiemann);
}

