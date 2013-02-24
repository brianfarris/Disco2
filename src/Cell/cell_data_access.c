#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdIO.h>
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
void cell_add_cons(struct Cell *oneCell, int q, double add){
  oneCell->cons[q] += add;
}
void cell_add_divB(struct Cell *oneCell, double add){
  oneCell->divB += add;
}
void cell_add_GradPsi(struct Cell *oneCell, int i, double add){
  oneCell->GradPsi[i] += add;
}
void cell_add_wiph(struct Cell *oneCell, double add){
  oneCell->wiph += add;
}
void cell_mult_psi(struct Cell *oneCell, double mult){
  oneCell->cons[PSI] *= mult;
}


void cell_set_prim(struct Cell ***theCells,int i,int j,int k,int q,double value) {
  theCells[k][i][j].prim[q] = value;
}
void cell_set_tiph(struct Cell ***theCells,int i,int j,int k,double value) {
  theCells[k][i][j].tiph = value;
}

void cell_printscreen(struct Cell ***theCells,struct Grid * theGrid){
  int i;

  for (i=0;i<(grid_Nghost_rmin(theGrid)+grid_Nghost_rmax(theGrid)+grid_N_r(theGrid));++i){
    double r = 0.5*(grid_r_faces(theGrid,i-1)+grid_r_faces(theGrid,i));
    printf("MyProc: %d r: %e i: %d rho: %e theCells[5][i][5].prim[PPP]: %e theCells[5][i][5].prim[URR]: %e theCells[5][i][5].prim[UPP]: %e theCells[5][i][5].prim[UZZ]: %e\n",MyProc,r,i,theCells[5][i][5].prim[RHO],theCells[5][i][5].prim[PPP],theCells[5][i][5].prim[URR],theCells[5][i][5].prim[UPP],theCells[5][i][5].prim[UZZ]);
  }

}
