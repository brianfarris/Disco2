#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Sim.h"
#include "../Headers/Face.h"
#include "../Headers/GravMass.h"
#include "../Headers/header.h"



void cell_set_wcell(struct Cell ***theCells,struct Sim *theSim){
  int i,j,k;
  for( k=0 ; k<sim_N_z(theSim) ; ++k ){
    for( i=0 ; i<sim_N_r(theSim) ; ++i ){
      double rp = sim_FacePos(theSim,i,R_DIR);
      double rm = sim_FacePos(theSim,i-1,R_DIR);
      double r = 0.5*(rm+rp);
      for( j=0 ; j<sim_N_p(theSim,i) ; ++j ){
        int jp = j+1;
        double w  = theCells[k][i][j ].prim[UPP];
        double wp = theCells[k][i][jp].prim[UPP];
        theCells[k][i][j].wiph = .5*r*(w+wp);
      }
    }
  }
}

void cell_set_wrigid(struct Cell ***theCells,struct Sim *theSim){
  int i,j,k;
  for( k=0 ; k<sim_N_z(theSim) ; ++k ){
    for( i=0 ; i<sim_N_r(theSim) ; ++i ){
      double w=0.0;
      double rp = sim_FacePos(theSim,i,R_DIR);
      double rm = sim_FacePos(theSim,i-1,R_DIR);
      double r = 0.5*(rm+rp);
      double Mring = 0.0;
      for( j=0 ; j<sim_N_p(theSim,i) ; ++j ){
        double vp = r*theCells[k][i][j].prim[UPP];
        double m = theCells[k][i][j].cons[DDD];
        Mring += m;
        w += vp*m;
      }
      for( j=0 ; j<sim_N_p(theSim,i) ; ++j ){
        theCells[k][i][j].wiph = w/Mring;       
      }
    }
  }

}


