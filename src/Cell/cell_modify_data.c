#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Sim.h"
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
//void cell_add_wiph(struct Cell *oneCell, double add){
//  oneCell->wiph += add;
//}
void cell_mult_psi(struct Cell *oneCell, double mult){
  oneCell->cons[PSI] *= mult;
}


void cell_set_prim(struct Cell ***theCells,int i,int j,int k,int q,double value) {
  theCells[k][i][j].prim[q] = value;
}
void cell_set_tiph(struct Cell ***theCells,int i,int j,int k,double value) {
  theCells[k][i][j].tiph = value;
}

