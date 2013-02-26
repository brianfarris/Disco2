#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Grid.h"
#include "../Headers/Face.h"
#include "../Headers/GravMass.h"
#include "../Headers/header.h"

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

/*
void cell_printscreen(struct Cell ***theCells,struct Grid * theGrid){
  int i;

  for (i=0;i<(grid_Nghost_rmin(theGrid)+grid_Nghost_rmax(theGrid)+grid_N_r(theGrid));++i){
    double r = 0.5*(grid_r_faces(theGrid,i-1)+grid_r_faces(theGrid,i));
    printf("MyProc: %d r: %e i: %d rho: %e theCells[5][i][5].prim[PPP]: %e theCells[5][i][5].prim[URR]: %e theCells[5][i][5].prim[UPP]: %e theCells[5][i][5].prim[UZZ]: %e\n",MyProc,r,i,theCells[5][i][5].prim[RHO],theCells[5][i][5].prim[PPP],theCells[5][i][5].prim[URR],theCells[5][i][5].prim[UPP],theCells[5][i][5].prim[UZZ]);
  }

}
*/
