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
    double zm = grid_z_faces(theGrid,k-1);
    double zp = grid_z_faces(theGrid,k);
    double dz = zp-zm;
    for( i=0 ; i<N_r_withghost ; ++i ){
      double rm = grid_r_faces(theGrid,i-1);
      double rp = grid_r_faces(theGrid,i);
      double dr = rp-rm;
      double r = .5*(rp+rm);
      for( j=0 ; j<grid_N_p(theGrid,i)-1 ; ++j ){
        struct Cell * cL = cell_single(theCells,i,j,k);
        struct Cell * cR = cell_single(theCells,i,j+1,k);
        riemann_driver( cL , cR ,theGrid, dr*dz , dt , r ,0.,0.,cell_dphi(cL),cell_dphi(cR),1);
      }
      struct Cell * cL = cell_single(theCells,i,grid_N_p(theGrid,i)-1,k);
      struct Cell * cR = cell_single(theCells,i,0,k);
      riemann_driver( cL , cR ,theGrid, dr*dz , dt , r,0.,0.,cell_dphi(cL),cell_dphi(cR),1 );
    }
  }
}


