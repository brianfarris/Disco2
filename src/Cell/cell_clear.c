#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdIO.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Grid.h"
#include "../Headers/Face.h"
#include "../Headers/GravMass.h"
#include "../Headers/header.h"


void cell_clear_w(struct Cell ***theCells,struct Grid * theGrid){
  int N_r_withghost = grid_N_r(theGrid)+grid_Nghost_rmin(theGrid)+grid_Nghost_rmax(theGrid);
  int N_z_withghost = grid_N_z(theGrid)+grid_Nghost_zmin(theGrid)+grid_Nghost_zmax(theGrid);
  int i,j,k,q;
  for( k=0 ; k<N_z_withghost ; ++k ){
    for( i=0 ; i<N_r_withghost ; ++i ){
      for( j=0 ; j<grid_N_p(theGrid,i) ; ++j ){
        theCells[k][i][j].wiph = 0.0;
      }
    }
  }
}

void cell_clear_divB( struct Cell *** theCells,struct Grid *theGrid ){
  int N_r_withghost = grid_N_r(theGrid)+grid_Nghost_rmin(theGrid)+grid_Nghost_rmax(theGrid);
  int N_z_withghost = grid_N_z(theGrid)+grid_Nghost_zmin(theGrid)+grid_Nghost_zmax(theGrid);
  int i,j,k;
  for( k=0 ; k<N_z_withghost ; ++k ){
    for( i=0 ; i<N_r_withghost ; ++i ){
      for( j=0 ; j<grid_N_p(theGrid,i) ; ++j ){
        theCells[k][i][j].divB = 0.0;
      }
    }
  }
}

void cell_clear_GradPsi( struct Cell *** theCells,struct Grid *theGrid){
  int N_r_withghost = grid_N_r(theGrid)+grid_Nghost_rmin(theGrid)+grid_Nghost_rmax(theGrid);
  int N_z_withghost = grid_N_z(theGrid)+grid_Nghost_zmin(theGrid)+grid_Nghost_zmax(theGrid);
  int i,j,k;
  for( k=0 ; k<N_z_withghost ; ++k ){
    for( i=0 ; i<N_r_withghost ; ++i ){
      for( j=0 ; j<grid_N_p(theGrid,i) ; ++j ){
        theCells[k][i][j].GradPsi[0] = 0.0;
        theCells[k][i][j].GradPsi[1] = 0.0;
        theCells[k][i][j].GradPsi[2] = 0.0;
      }
    }
  }
}



