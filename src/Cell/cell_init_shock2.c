#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Sim.h"
#include "../Headers/Face.h"
#include "../Headers/GravMass.h"
#include "../Headers/header.h"

void cell_single_init_shock2(struct Cell *theCell, struct Sim *theSim,int i,int j,int k){
  double rhoL  = 1.0;
  double rhoR  = 0.125;
  double vL = 0.0;
  double vR = 0.0;
  double pL = 1.0;
  double pR = 0.1;

  double rm = sim_FacePos(theSim,i-1,R_DIR);
  double rp = sim_FacePos(theSim,i,R_DIR);
  double r = 0.5*(rm+rp);
  double zm = sim_FacePos(theSim,k-1,Z_DIR);
  double zp = sim_FacePos(theSim,k,Z_DIR);
  double z = 0.5*(zm+zp);
  double t = theCell->tiph-.5*theCell->dphi;

  if(r*cos(t) < 0.0)
  {
    theCell->prim[RHO] = rhoL;
    theCell->prim[PPP] = pL;
    theCell->prim[URR] = vL * cos(t);
    theCell->prim[UPP] = -vL * sin(t) / r;
    theCell->prim[UZZ] = 0.0;
    theCell->divB = 0.0;
    theCell->GradPsi[0] = 0.0;
    theCell->GradPsi[1] = 0.0;
    theCell->GradPsi[2] = 0.0;
  }
  else
  {
    theCell->prim[RHO] = rhoR;
    theCell->prim[PPP] = pR;
    theCell->prim[URR] = vL * cos(t);
    theCell->prim[UPP] = -vL * sin(t) / r;
    theCell->prim[UZZ] = 0.0;
    theCell->divB = 0.0;
    theCell->GradPsi[0] = 0.0;
    theCell->GradPsi[1] = 0.0;
    theCell->GradPsi[2] = 0.0;
  }

  //TODO: Not sure what this is for.  Ask someone if important.
  //if(sim_NUM_C(theSim)<sim_NUM_Q(theSim)) theCell->prim[sim_NUM_C(theSim)] = Qq;
}

void cell_init_shock2(struct Cell ***theCells,struct Sim *theSim,struct MPIsetup * theMPIsetup) {

  double rhoL  = 1.0;
  double rhoR  = 0.125;
  double vL = 0.0;
  double vR = 0.0;
  double pL = 1.0;
  double pR = 0.1;

  int i, j, k;
  for (k = 0; k < sim_N(theSim,Z_DIR); k++) {
    for (i = 0; i < sim_N(theSim,R_DIR); i++) {
      double rm = sim_FacePos(theSim,i-1,R_DIR);
      double rp = sim_FacePos(theSim,i,R_DIR);
      double r = 0.5*(rm+rp);
      double zm = sim_FacePos(theSim,k-1,Z_DIR);
      double zp = sim_FacePos(theSim,k,Z_DIR);
      double z = 0.5*(zm+zp);
      for (j = 0; j < sim_N_p(theSim,i); j++) {
        double t = theCells[k][i][j].tiph-.5*theCells[k][i][j].dphi;
        if(r*cos(t) < 0.0) {
          theCells[k][i][j].prim[RHO] = rhoL;
          theCells[k][i][j].prim[PPP] = pL;
          theCells[k][i][j].prim[URR] = vL * cos(t);
          theCells[k][i][j].prim[UPP] = -vL * sin(t) / r;
          theCells[k][i][j].prim[UZZ] = 0.0;
          theCells[k][i][j].divB = 0.0;
          theCells[k][i][j].GradPsi[0] = 0.0;
          theCells[k][i][j].GradPsi[1] = 0.0;
          theCells[k][i][j].GradPsi[2] = 0.0;
        } else {
          theCells[k][i][j].prim[RHO] = rhoR;
          theCells[k][i][j].prim[PPP] = pR;
          theCells[k][i][j].prim[URR] = vL * cos(t);
          theCells[k][i][j].prim[UPP] = -vL * sin(t) / r;
          theCells[k][i][j].prim[UZZ] = 0.0;
          theCells[k][i][j].divB = 0.0;
          theCells[k][i][j].GradPsi[0] = 0.0;
          theCells[k][i][j].GradPsi[1] = 0.0;
          theCells[k][i][j].GradPsi[2] = 0.0;
        }
      }
    }
  }
}
