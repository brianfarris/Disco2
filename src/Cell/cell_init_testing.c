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
  double rm = sim_FacePos(theSim,i-1,R_DIR);
  double rp = sim_FacePos(theSim,i,R_DIR);
  double r = 0.5*(rm+rp);


  double phi = theCell->tiph - 0.5*theCell->dphi;
  double rho = 1.;///(1.+r); 
  double omega = 100.0; 
  double vr = 0.0;//r;
  if( rho<sim_RHO_FLOOR(theSim) ) rho = sim_RHO_FLOOR(theSim);

  double Pp = rho/sim_GAMMALAW(theSim);

  theCell->prim[RHO] = rho;
  theCell->prim[PPP] = Pp;
  theCell->prim[URR] = vr;
  theCell->prim[UPP] = omega;
  theCell->prim[UZZ] = 0.0;
  theCell->wiph = 0.0;
  theCell->divB = 0.0;
  theCell->GradPsi[0] = 0.0;
  theCell->GradPsi[1] = 0.0;
  theCell->GradPsi[2] = 0.0;
  theCell->prim[BRR] = 0.0;
  theCell->prim[BPP] = 0.0;
  theCell->prim[BZZ] = 0.0;
  theCell->prim[PSI] = 0.0;
  theCell->prim[ARR] = 0.0;
  theCell->prim[APP] = 0.0;
  theCell->prim[AZZ] = 0.0;
  theCell->prim[PHI] = 0.0;
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
        double rho = 1.;///(1.+r); 
        double omega = 100.0; 
        double vr = 0.0;//r;
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
        theCells[k][i][j].prim[BRR] = 0.0;
        theCells[k][i][j].prim[BPP] = 0.0;
        theCells[k][i][j].prim[BZZ] = 0.0;
        theCells[k][i][j].prim[PSI] = 0.0;
        theCells[k][i][j].prim[ARR] = 0.0;
        theCells[k][i][j].prim[APP] = 0.0;
        theCells[k][i][j].prim[AZZ] = 0.0;
        theCells[k][i][j].prim[PHI] = 0.0;
      }
    }
  }

}
