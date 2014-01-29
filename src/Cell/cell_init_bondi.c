#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Sim.h"
#include "../Headers/Face.h"
#include "../Headers/GravMass.h"
#include "../Headers/MPIsetup.h"
#include "../Headers/header.h"

void cell_single_init_bondi(struct Cell *theCell, struct Sim *theSim,int i,int j,int k){
  double GAMMALAW = sim_GAMMALAW(theSim);
  double RG = sim_GravRadius(theSim);
  double rho = sim_InitPar1(theSim);
  double P = sim_InitPar2(theSim);

  theCell->prim[RHO] = rho;
  theCell->prim[PPP] = P;
  theCell->prim[URR] = 0.0;
  theCell->prim[UPP] = 0.0; 
  theCell->prim[UZZ] = 0.0;
  theCell->wiph = 0.0;

}

void cell_init_bondi(struct Cell ***theCells,struct Sim *theSim,struct MPIsetup * theMPIsetup) {
  double GAMMALAW = sim_GAMMALAW(theSim);
  double RG = sim_GravRadius(theSim);
  double rho = sim_InitPar1(theSim);
  double P = sim_InitPar2(theSim);
  
  int i, j,k;
  for (k = 0; k < sim_N(theSim,Z_DIR); k++) {
    for (i = 0; i < sim_N(theSim,R_DIR); i++) {
      for (j = 0; j < sim_N_p(theSim,i); j++) {
        theCells[k][i][j].prim[RHO] = rho;
        theCells[k][i][j].prim[PPP] = P;
        theCells[k][i][j].prim[URR] = 0.0;
        theCells[k][i][j].prim[UPP] = 0.0; 
        theCells[k][i][j].prim[UZZ] = 0.0;
        theCells[k][i][j].wiph = 0.0;
      }
    }
  }

}
