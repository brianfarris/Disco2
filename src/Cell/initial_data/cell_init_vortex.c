#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../../Headers/Cell.h"
#include "../../Headers/Sim.h"
#include "../../Headers/Face.h"
#include "../../Headers/GravMass.h"
#include "../../Headers/header.h"

void cell_single_init_vortex(struct Cell *theCell, struct Sim *theSim,int i,int j,int k){
  double rho  = 1.0;
  double P0   = 1.0;
  double R    = 1.0;
  double A    = 0.5;
  double g    = 3.0;

  double rm = sim_FacePos(theSim,i-1,R_DIR);
  double rp = sim_FacePos(theSim,i,R_DIR);
  double r = 0.5*(rm+rp);
  double Pp,vp,Qq;
  double t = theCell->tiph-.5*theCell->dphi;
  if( r >= R ){
    Pp = P0;
    vp = 0.0;
  }else{
    Pp = P0*( 1. - A*exp( g*( r/R - R/(R-r) + 1. ) ) );
    vp = sqrt( A*g*P0/rho*(2.*R-r)/R*r*r/pow(R-r,2.)*exp( g*( r/R - R/(R-r) + 1. ) ) );
    if( cos(t) < 0.0 ) Qq = 0.0; else Qq = 1.0;
    theCell->prim[RHO] = rho;
    theCell->prim[PPP] = Pp;
    theCell->prim[URR] = 0.0;
    theCell->prim[UPP] = vp/r;
    theCell->prim[UZZ] = 0.0;
    theCell->divB = 0.0;
    theCell->GradPsi[0] = 0.0;
    theCell->GradPsi[1] = 0.0;
    theCell->GradPsi[2] = 0.0;
    if(sim_NUM_C(theSim)<sim_NUM_Q(theSim)) theCell->prim[sim_NUM_C(theSim)] = Qq;
  }
}

void cell_init_vortex(struct Cell ***theCells,struct Sim *theSim,struct MPIsetup * theMPIsetup) {

  double rho  = 1.0;
  double P0   = 1.0;
  double R    = 1.0;
  double A    = 0.5;
  double g    = 3.0;

  int i, j,k;
  for (k = 0; k < sim_N(theSim,Z_DIR); k++) {
    for (i = 0; i < sim_N(theSim,R_DIR); i++) {
      double rm = sim_FacePos(theSim,i-1,R_DIR);
      double rp = sim_FacePos(theSim,i,R_DIR);
      double r = 0.5*(rm+rp);
      for (j = 0; j < sim_N_p(theSim,i); j++) {
        double Pp,vp,Qq;
        double t = theCells[k][i][j].tiph-.5*theCells[k][i][j].dphi;
        if( r >= R ){
          Pp = P0;
          vp = 0.0;
        }else{
          Pp = P0*( 1. - A*exp( g*( r/R - R/(R-r) + 1. ) ) );
          vp = sqrt( A*g*P0/rho*(2.*R-r)/R*r*r/pow(R-r,2.)*exp( g*( r/R - R/(R-r) + 1. ) ) );
          if( cos(t) < 0.0 ) Qq = 0.0; else Qq = 1.0;
          theCells[k][i][j].prim[RHO] = rho;
          theCells[k][i][j].prim[PPP] = Pp;
          theCells[k][i][j].prim[URR] = 0.0;
          theCells[k][i][j].prim[UPP] = vp/r;
          theCells[k][i][j].prim[UZZ] = 0.0;
          theCells[k][i][j].divB = 0.0;
          theCells[k][i][j].GradPsi[0] = 0.0;
          theCells[k][i][j].GradPsi[1] = 0.0;
          theCells[k][i][j].GradPsi[2] = 0.0;
          if(sim_NUM_C(theSim)<sim_NUM_Q(theSim)) theCells[k][i][j].prim[sim_NUM_C(theSim)] = Qq;

        }
      }
    }
  }
}
