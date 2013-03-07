#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Grid.h"
#include "../Headers/Face.h"
#include "../Headers/GravMass.h"
#include "../Headers/header.h"


void cell_clear_w(struct Cell ***theCells,struct Grid * theGrid){
  int i,j,k,q;
  for( k=0 ; k<grid_N_z(theGrid) ; ++k ){
    for( i=0 ; i<grid_N_r(theGrid) ; ++i ){
      for( j=0 ; j<grid_N_p(theGrid,i) ; ++j ){
        theCells[k][i][j].wiph = 0.0;
      }
    }
  }
}

void cell_clear_divB( struct Cell *** theCells,struct Grid *theGrid ){
  int i,j,k;
  for( k=0 ; k<grid_N_z(theGrid) ; ++k ){
    for( i=0 ; i<grid_N_r(theGrid) ; ++i ){
      for( j=0 ; j<grid_N_p(theGrid,i) ; ++j ){
        theCells[k][i][j].divB = 0.0;
      }
    }
  }
}

void cell_clear_GradPsi( struct Cell *** theCells,struct Grid *theGrid){
  int i,j,k;
  for( k=0 ; k<grid_N_z(theGrid) ; ++k ){
    for( i=0 ; i<grid_N_r(theGrid) ; ++i ){
      for( j=0 ; j<grid_N_p(theGrid,i) ; ++j ){
        theCells[k][i][j].GradPsi[0] = 0.0;
        theCells[k][i][j].GradPsi[1] = 0.0;
        theCells[k][i][j].GradPsi[2] = 0.0;
      }
    }
  }
}



