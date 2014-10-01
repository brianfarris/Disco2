#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Sim.h"
#include "../Headers/Face.h"
#include "../Headers/GravMass.h"
#include "../Headers/header.h"

double cell_prim(struct Cell *theCell,int q){
  return(theCell->prim[q]);
}
double cell_grad(struct Cell *theCell,int q){
  return(theCell->grad[q]);
}
double cell_gradp(struct Cell *theCell,int q){
  return(theCell->gradp[q]);
}
struct Cell *cell_single(struct Cell ***theCells,int i,int j,int k){
  return &(theCells[k][i][j]);
}
double cell_tiph(struct Cell *oneCell){
  return(oneCell->tiph);
}
double cell_dphi(struct Cell *oneCell){
  return(oneCell->dphi);
}
double cell_wiph(struct Cell *oneCell){
  return(oneCell->wiph);
}
double cell_Cool(struct Cell ***theCells,int i,int j,int k){
  return(theCells[k][i][j].Cool);
}
