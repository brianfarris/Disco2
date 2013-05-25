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
  theRiemann->UL = malloc(sizeof(double)*sim_NUM_Q(theSim));
  theRiemann->UR = malloc(sizeof(double)*sim_NUM_Q(theSim));
  theRiemann->Ustar = malloc(sizeof(double)*sim_NUM_Q(theSim));
  theRiemann->FL = malloc(sizeof(double)*sim_NUM_Q(theSim));
  theRiemann->FR = malloc(sizeof(double)*sim_NUM_Q(theSim));
  theRiemann->Fstar = malloc(sizeof(double)*sim_NUM_Q(theSim));
  theRiemann->F = malloc(sizeof(double)*sim_NUM_Q(theSim));
  int q;
  for (q=0;q<sim_NUM_Q(theSim);++q){
    theRiemann->primL[q]=0.;
    theRiemann->primR[q]=0.;
    theRiemann->UL[q]=0.;
    theRiemann->UR[q]=0.;
    theRiemann->Ustar[q]=0.;
    theRiemann->FL[q]=0.;
    theRiemann->FR[q]=0.;
    theRiemann->Fstar[q]=0.;
    theRiemann->F[q]=0.;
  }
  int iter;
  for (iter=0;iter<3;++iter){
    theRiemann->n[iter]=0; //initialize
  }
  return(theRiemann);
}

void riemann_destroy(struct Riemann * theRiemann){
  free(theRiemann->F);
  free(theRiemann->Fstar);
  free(theRiemann->FR);
  free(theRiemann->FL);
  free(theRiemann->primL);
  free(theRiemann->primR);
  free(theRiemann->UR);
  free(theRiemann->UL);
  free(theRiemann->Ustar);
  free(theRiemann);
}

