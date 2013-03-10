#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Sim.h"
#include "../Headers/Face.h"
#include "../Headers/GravMass.h"
#include "../Headers/header.h"


void cell_clear_w(struct Cell ***theCells,struct Sim * theSim){
  int i,j,k,q;
  for( k=0 ; k<sim_N_z(theSim) ; ++k ){
    for( i=0 ; i<sim_N_r(theSim) ; ++i ){
      for( j=0 ; j<sim_N_p(theSim,i) ; ++j ){
        theCells[k][i][j].wiph = 0.0;
      }
    }
  }
}

void cell_clear_divB( struct Cell *** theCells,struct Sim *theSim ){
  int i,j,k;
  for( k=0 ; k<sim_N_z(theSim) ; ++k ){
    for( i=0 ; i<sim_N_r(theSim) ; ++i ){
      for( j=0 ; j<sim_N_p(theSim,i) ; ++j ){
        theCells[k][i][j].divB = 0.0;
      }
    }
  }
}

void cell_clear_GradPsi( struct Cell *** theCells,struct Sim *theSim){
  int i,j,k;
  for( k=0 ; k<sim_N_z(theSim) ; ++k ){
    for( i=0 ; i<sim_N_r(theSim) ; ++i ){
      for( j=0 ; j<sim_N_p(theSim,i) ; ++j ){
        theCells[k][i][j].GradPsi[0] = 0.0;
        theCells[k][i][j].GradPsi[1] = 0.0;
        theCells[k][i][j].GradPsi[2] = 0.0;
      }
    }
  }
}



