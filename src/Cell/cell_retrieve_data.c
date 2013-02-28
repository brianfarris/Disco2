#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Grid.h"
#include "../Headers/Face.h"
#include "../Headers/GravMass.h"
#include "../Headers/header.h"

double *cell_get_prims(struct Cell *theCell){
  return(theCell->prim);
}
double *cell_single_get_grad(struct Cell *theCell){
  return(theCell->grad);
}
double *cell_single_get_gradp(struct Cell *theCell){
  return(theCell->gradp);
}
struct Cell *cell_pointer(struct Cell ***theCells,int i,int j,int k){
  return &(theCells[k][i][j]);
}
double cell_single_tiph(struct Cell *oneCell){
  return(oneCell->tiph);
}
double cell_single_dphi(struct Cell *oneCell){
  return(oneCell->dphi);
}
double cell_single_wiph(struct Cell *oneCell){
  return(oneCell->wiph);
}

