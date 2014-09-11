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

void cell_single_init_mhdexp(struct Cell *theCell, struct Sim *theSim,int i,int j,int k){
  printf("ERROR, DO NOT CALL THIS FOR MHDEXP\n");
  exit(0);
}

void cell_init_mhdexp(struct Cell ***theCells,struct Sim *theSim,struct MPIsetup * theMPIsetup) {
  double GAMMALAW = sim_GAMMALAW(theSim);

  double R0 = 0.1;
  double B0 = 1.0;

  int i, j,k;
  for (k = 0; k < sim_N(theSim,Z_DIR); k++) {
    for (i = 0; i < sim_N(theSim,R_DIR); i++) {
      double rm = sim_FacePos(theSim,i-1,R_DIR);
      double rp = sim_FacePos(theSim,i,R_DIR);
      double r = 0.5*(rm+rp);
      for (j = 0; j < sim_N_p(theSim,i); j++) {
        double phi = cell_tiph(cell_single(theCells,i,j,k)) - 0.5*cell_dphi(cell_single(theCells,i,j,k))-M_PI/4.;
        theCells[k][i][j].prim[RHO] = 1.0;
        theCells[k][i][j].prim[PPP] = 0.1;
        theCells[k][i][j].prim[URR] = 0.0;
        theCells[k][i][j].prim[UPP] = 0.0;
        theCells[k][i][j].prim[UZZ] = 0.0;
        theCells[k][i][j].prim[BRR] = B0*cos(phi);
        theCells[k][i][j].prim[BPP] = -B0*sin(phi);
        theCells[k][i][j].prim[BZZ] = 0.0;
        theCells[k][i][j].prim[PSI] = 0.0;
        theCells[k][i][j].wiph = 0.0;
        theCells[k][i][j].divB = 0.0;
        theCells[k][i][j].GradPsi[0] = 0.0;
        theCells[k][i][j].GradPsi[1] = 0.0;
        theCells[k][i][j].GradPsi[2] = 0.0;
        if (r<R0){
          theCells[k][i][j].prim[PPP] = 10.0;
        }
      }
    }
  }

}

