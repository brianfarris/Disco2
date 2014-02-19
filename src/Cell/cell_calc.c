#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Sim.h"
#include "../Headers/Face.h"
#include "../Headers/GravMass.h"
#include "../Headers/header.h"

void cell_calc_cons( struct Cell *** theCells,struct Sim *theSim ){
  double GAMMALAW = sim_GAMMALAW(theSim);

  int i,j,k;
  double pos[3];
  for( k=0 ; k<sim_N(theSim,Z_DIR) ; ++k ){
    double zp = sim_FacePos(theSim,k,Z_DIR);
    double zm = sim_FacePos(theSim,k-1,Z_DIR);
    double dz = zp-zm;
    pos[Z_DIR] = 0.5*(zm+zp);
    for( i=0 ; i<sim_N(theSim,R_DIR) ; ++i ){
      double rp = sim_FacePos(theSim,i,R_DIR);
      double rm = sim_FacePos(theSim,i-1,R_DIR);
      double r = .5*(rp+rm);
      pos[R_DIR] = 0.5*(rm+rp);
      for( j=0 ; j<sim_N_p(theSim,i) ; ++j ){
        struct Cell * c = &(theCells[k][i][j]);
        double dV = .5*(rp*rp - rm*rm)*c->dphi*dz;
        pos[P_DIR] = c->tiph - 0.5*c->dphi;
        cell_prim2cons( c->prim , c->cons , pos , dV,theSim);
      }    
    }    
  }
}

void cell_calc_prim( struct Cell *** theCells ,struct Sim * theSim){
  int i,j,k;
  double pos[3];
  for( k=0 ; k<sim_N(theSim,Z_DIR) ; ++k ){
    double zm = sim_FacePos(theSim,k-1,Z_DIR);
    double zp = sim_FacePos(theSim,k,Z_DIR);
    double dz = zp-zm;
    pos[Z_DIR] = 0.5*(zm+zp);
    for( i=0 ; i<sim_N(theSim,R_DIR) ; ++i ){
      double rm = sim_FacePos(theSim,i-1,R_DIR);
      double rp = sim_FacePos(theSim,i,R_DIR);
      pos[R_DIR] = 0.5*(rm+rp);
      for( j=0 ; j<sim_N_p(theSim,i) ; ++j ){
        struct Cell * c = cell_single(theCells,i,j,k);
        double dV = .5*(rp*rp-rm*rm)*c->dphi*dz;
        pos[P_DIR] = c->tiph - 0.5*c->dphi;
        cell_cons2prim( c->cons , c->prim , pos , dV ,theSim);
      }
    }
  }
}

