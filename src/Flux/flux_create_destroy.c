#define FLUX_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Flux.h"
#include "../Headers/header.h"

struct Flux *flux_create(struct Grid *theGrid){
  struct Flux * theFlux = malloc(sizeof(struct Flux));
  Flux->primL = malloc(sizeof(double)*grid_NUM_Q(theGrid));
  Flux->primR = malloc(sizeof(double)*grid_NUM_Q(theGrid));
  Flux->Uk = malloc(sizeof(double)*grid_NUM_Q(theGrid));
  Flux->Ustar = malloc(sizeof(double)*grid_NUM_Q(theGrid));
  Flux->F = malloc(sizeof(double)*grid_NUM_Q(theGrid));
  int q;
  for (q=0;q<NUM_Q;++q){
    Flux->primL[q]=0.;
    Flux->primR[q]=0.;
    Flux->Uk[q]=0.;
    Flux->Ustar[q]=0.;
     Flux->F[q]=0.;
  }
  return(theFlux);
}

void flux_destroy(struct Flux * theFlux){
  free(theFlux->F);
  free(theFlux->primL);
  free(theFlux->primR);
  free(theFlux->Uk);
  free(theFlux->Ustar);
  free(theFlux);
}

