#define RIEMANN_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Grid.h"
#include "../Headers/Riemann.h"
#include "../Headers/header.h"

struct Riemann *riemann_create(struct Grid *theGrid){
  struct Riemann * theRiemann = malloc(sizeof(struct Riemann));
  theRiemann->primL = malloc(sizeof(double)*grid_NUM_Q(theGrid));
  theRiemann->primR = malloc(sizeof(double)*grid_NUM_Q(theGrid));
  theRiemann->Uk = malloc(sizeof(double)*grid_NUM_Q(theGrid));
  theRiemann->Ustar = malloc(sizeof(double)*grid_NUM_Q(theGrid));
  theRiemann->F = malloc(sizeof(double)*grid_NUM_Q(theGrid));
  int q;
  for (q=0;q<grid_NUM_Q(theGrid);++q){
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

