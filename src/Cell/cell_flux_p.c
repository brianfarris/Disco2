#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Riemann.h"
#include "../Headers/Grid.h"
#include "../Headers/Face.h"
#include "../Headers/GravMass.h"
#include "../Headers/header.h"

void cell_flux_p( struct Cell *** theCells ,struct Grid *theGrid, double dt ){
  int N_r_withghost = grid_N_r(theGrid)+grid_Nghost_rmin(theGrid)+grid_Nghost_rmax(theGrid);
  int N_z_withghost = grid_N_z(theGrid)+grid_Nghost_zmin(theGrid)+grid_Nghost_zmax(theGrid);
  int i,j,k;
  for( k=0 ; k<N_z_withghost ; ++k ){
    for( i=0 ; i<N_r_withghost ; ++i ){
      for( j=0 ; j<grid_N_p(theGrid,i)-1 ; ++j ){
        struct Riemann * theRiemann = riemann_create(theGrid);
        riemann_setup_p(theRiemann,theCells,theGrid,i,j,j+1,k);
        riemann_hllc( theRiemann,theGrid,dt,1);
        riemann_destroy(theRiemann);
      }
      struct Riemann * theRiemann = riemann_create(theGrid);
      riemann_setup_p(theRiemann,theCells,theGrid,i,grid_N_p(theGrid,i)-1,0,k);
      riemann_hllc(theRiemann, theGrid,dt,1 );
      riemann_destroy(theRiemann);
    }
  }
}


