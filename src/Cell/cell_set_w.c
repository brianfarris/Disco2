#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Grid.h"
#include "../Headers/Face.h"
#include "../Headers/GravMass.h"
#include "../Headers/header.h"



void cell_set_wcell(struct Cell ***theCells,struct Grid *theGrid){
  int N_r_withghost = grid_N_r(theGrid)+grid_Nghost_rmin(theGrid)+grid_Nghost_rmax(theGrid);
  int N_z_withghost = grid_N_z(theGrid)+grid_Nghost_zmin(theGrid)+grid_Nghost_zmax(theGrid);
  int i,j,k;
  for( k=0 ; k<N_z_withghost ; ++k ){
    for( i=0 ; i<N_r_withghost ; ++i ){
      double rp = grid_r_faces(theGrid,i);
      double rm = grid_r_faces(theGrid,i-1);
      double r = 0.5*(rm+rp);
      for( j=0 ; j<grid_N_p(theGrid,i) ; ++j ){
        int jp = j+1;
        double w  = theCells[k][i][j ].prim[UPP];
        double wp = theCells[k][i][jp].prim[UPP];
        theCells[k][i][j].wiph = .5*r*(w+wp);
      }
    }
  }
}

void cell_set_wrigid(struct Cell ***theCells,struct Grid *theGrid){
  int N_r_withghost = grid_N_r(theGrid)+grid_Nghost_rmin(theGrid)+grid_Nghost_rmax(theGrid);
  int N_z_withghost = grid_N_z(theGrid)+grid_Nghost_zmin(theGrid)+grid_Nghost_zmax(theGrid);
  int i,j,k;
  for( k=0 ; k<N_z_withghost ; ++k ){
    for( i=0 ; i<N_r_withghost ; ++i ){
      double w=0.0;
      double rp = grid_r_faces(theGrid,i);
      double rm = grid_r_faces(theGrid,i-1);
      double r = 0.5*(rm+rp);
      double Mring = 0.0;
      for( j=0 ; j<grid_N_p(theGrid,i) ; ++j ){
        double vp = r*theCells[k][i][j].prim[UPP];
        double m = theCells[k][i][j].cons[DDD];
        Mring += m;
        w += vp*m;
      }
      for( j=0 ; j<grid_N_p(theGrid,i) ; ++j ){
        theCells[k][i][j].wiph = w/Mring;       
        //printf("theCells[%d][%d][%d].wiph: %e\n",k,i,j,theCells[k][i][j].wiph);
      }
    }
  }

}


