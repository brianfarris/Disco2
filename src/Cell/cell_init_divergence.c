#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Sim.h"
#include "../Headers/Face.h"
#include "../Headers/GravMass.h"
#include "../Headers/header.h"

void cell_single_init_divergence(struct Cell *theCell, struct Sim *theSim,int i,int j,int k){

  double B0 = 0.01;
  double Rho0 = 1.0;
  double P0 = 1.0;
  double xinner = 1.0;
  double x0 = 0.3;

  double rm = sim_FacePos(theSim,i-1,R_DIR);
  double rp = sim_FacePos(theSim,i,R_DIR);
  double r = 0.5*(rm+rp);
  double Pp = P0;

  double phi = theCell->tiph - .5*theCell->dphi;

  double x = r*cos(phi)-x0;
  double y = r*sin(phi);

  double rl = sqrt(x*x+y*y);

  double Br;
  if (rl<0.15){
    Br = B0*rl;
  } else{
    Br = 0.0;
  }

  double Bx = Br*x/rl;
  double By = Br*y/rl;

  theCell->prim[RHO] = Rho0;
  theCell->prim[PPP] = Pp;
  theCell->prim[URR] = 0.0;
  theCell->prim[UPP] = 0.0;
  theCell->prim[UZZ] = 0.0;
  theCell->prim[BRR] =  Bx*cos(phi) + By*sin(phi);
  theCell->prim[BPP] = -Bx*sin(phi) + By*cos(phi);
  theCell->prim[BZZ] = 0.0;
  theCell->prim[PSI] = 0.0;

  theCell->divB = 0.0;
  theCell->GradPsi[0] = 0.0;
  theCell->GradPsi[1] = 0.0;
  theCell->GradPsi[2] = 0.0;
}

void cell_init_divergence(struct Cell ***theCells,struct Sim *theSim,struct MPIsetup * theMPIsetup) {

  double B0 = 0.01;
  double Rho0 = 1.0;
  double P0 = 1.0;
  double xinner = 1.0;
  double x0 = 0.3;

  int i, j,k;
  for (k = 0; k < sim_N(theSim,Z_DIR); k++) {
    for (i = 0; i < sim_N(theSim,R_DIR); i++) {
      double rm = sim_FacePos(theSim,i-1,R_DIR);
      double rp = sim_FacePos(theSim,i,R_DIR);
      double r = 0.5*(rm+rp);
      double Pp = P0;

      for (j = 0; j < sim_N_p(theSim,i); j++) {
        double phi = theCells[k][i][j].tiph - .5*theCells[k][i][j].dphi;

        double x = r*cos(phi)-x0;
        double y = r*sin(phi);

        double rl = sqrt(x*x+y*y);

        double Br; 
        if (rl<0.15){
          Br = B0*rl;
        } else{
          Br = 0.0;
        }

        double Bx = Br*x/rl;
        double By = Br*y/rl;

        theCells[k][i][j].prim[RHO] = Rho0;
        theCells[k][i][j].prim[PPP] = Pp;
        theCells[k][i][j].prim[URR] = 0.0;
        theCells[k][i][j].prim[UPP] = 0.0;
        theCells[k][i][j].prim[UZZ] = 0.0;
        theCells[k][i][j].prim[BRR] =  Bx*cos(phi) + By*sin(phi);
        theCells[k][i][j].prim[BPP] = -Bx*sin(phi) + By*cos(phi);
        theCells[k][i][j].prim[BZZ] = 0.0;
        theCells[k][i][j].prim[PSI] = 0.0;

        theCells[k][i][j].divB = 0.0;
        theCells[k][i][j].GradPsi[0] = 0.0;
        theCells[k][i][j].GradPsi[1] = 0.0;
        theCells[k][i][j].GradPsi[2] = 0.0;
      }
    }
  }
}
