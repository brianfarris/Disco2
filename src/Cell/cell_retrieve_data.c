#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Grid.h"
#include "../Headers/Face.h"
#include "../Headers/GravMass.h"
#include "../Headers/header.h"

 double *cell_get_prims(struct Cell ***theCells,int i,int j,int k){
  return(theCells[k][i][j].prim);
}
double *cell_get_cons(struct Cell ***theCells,int i,int j,int k){
  return(theCells[k][i][j].cons);
}

double *cell_single_get_prims(struct Cell *theCells){
  return(theCells->prim);
}
double *cell_single_get_grad(struct Cell *theCells){
  return(theCells->grad);
}
double *cell_single_get_gradp(struct Cell *theCells){
  return(theCells->gradp);
}

struct Cell *cell_pointer(struct Cell ***theCells,int i,int j,int k){
  return &(theCells[k][i][j]);
}

double cell_tiph(struct Cell ***theCells,int i,int j,int k){
  return(theCells[k][i][j].tiph);
}
double cell_single_tiph(struct Cell *oneCell){
  return(oneCell->tiph);
}

double cell_dphi(struct Cell ***theCells,int i,int j,int k){
  return(theCells[k][i][j].dphi);
}
double cell_single_dphi(struct Cell *oneCell){
  return(oneCell->dphi);
}
double cell_single_wiph(struct Cell *oneCell){
  return(oneCell->wiph);
}

