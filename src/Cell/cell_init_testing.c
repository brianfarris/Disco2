#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Sim.h"
#include "../Headers/Face.h"
#include "../Headers/GravMass.h"
#include "../Headers/header.h"

void cell_single_init_testing(struct Cell *theCell, struct Sim *theSim,int i,int j,int k){
  printf("ERROR. cell_single_init_shear isnt set up right now\n");
  exit(0);
}

void cell_init_testing(struct Cell ***theCells,struct Sim *theSim,struct MPIsetup * theMPIsetup) {
  int i, j,k;
  for (k = 0; k < sim_N(theSim,Z_DIR); k++) {
    for (i = 0; i < sim_N(theSim,R_DIR); i++) {
      double rm = sim_FacePos(theSim,i-1,R_DIR);
      double rp = sim_FacePos(theSim,i,R_DIR);
      double r = 0.5*(rm+rp);

      for (j = 0; j < sim_N_p(theSim,i); j++) {

        double phi = theCells[k][i][j].tiph - 0.5*theCells[k][i][j].dphi;
        double rho = 1./(1.+r); 
        double omega = 0.0; 
        double vr = r;
        if( rho<sim_RHO_FLOOR(theSim) ) rho = sim_RHO_FLOOR(theSim);
  
        double Pp = rho/sim_GAMMALAW(theSim);

        theCells[k][i][j].prim[RHO] = rho;
        theCells[k][i][j].prim[PPP] = Pp;
        theCells[k][i][j].prim[URR] = vr;
        theCells[k][i][j].prim[UPP] = omega;
        theCells[k][i][j].prim[UZZ] = 0.0;
        theCells[k][i][j].wiph = 0.0;
        theCells[k][i][j].divB = 0.0;
        theCells[k][i][j].GradPsi[0] = 0.0;
        theCells[k][i][j].GradPsi[1] = 0.0;
        theCells[k][i][j].GradPsi[2] = 0.0;
      }
    }
  }

}
