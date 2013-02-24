#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../Headers/cell.h"
#include "../Headers/grid.h"
#include "../Headers/face.h"
#include "../Headers/GravMass.h"
#include "../Headers/onestep.h"
#include "../Headers/header.h"

void onestep(struct Cell *** theCells,struct Grid * theGrid,struct GravMass * theGravMasses,double dt,double RK){
  int N_r_withghost = grid_N_r(theGrid)+grid_Nghost_rmin(theGrid)+grid_Nghost_rmax(theGrid);
  int N_z_withghost = grid_N_z(theGrid)+grid_Nghost_zmin(theGrid)+grid_Nghost_zmax(theGrid);
  int *nri = malloc(N_r_withghost*sizeof(int));
  int *nzk = malloc(N_z_withghost*sizeof(int));
   //onestep
  int Nfr,Nfz;
  struct Face *theFaces_r = face_create_r(theCells,theGrid,&Nfr,nri);
  struct Face *theFaces_z = face_create_z(theCells,theGrid,&Nfz,nzk);
  cell_clean_pi(theCells,theGrid);
  cell_clear_w(theCells,theGrid);
  if( MOVE_CELLS == C_WCELL ) cell_set_wcell( theCells ,theGrid);
  if( MOVE_CELLS == C_RIGID ) cell_set_wrigid( theCells ,theGrid);
  cell_adjust_RK_cons( theCells, theGrid, RK);

  cell_clear_divB(theCells,theGrid);
  cell_clear_GradPsi(theCells,theGrid);
  //Phi Flux
  cell_plm_p( theCells,theGrid );
  cell_flux_p( theCells ,theGrid, dt );
  //R Flux
  cell_plm_rz( theCells ,theGrid, theFaces_r , Nfr , 0 );
  int n;
  for( n=0 ; n<Nfr ; ++n ){
    face_riemann_r( face_pointer(theFaces_r,n) , dt );
  }

  //Z Flux
  if( N_z_global != 1 ){
    cell_plm_rz( theCells ,theGrid, theFaces_z , Nfz , 1 );
    for( n=0 ; n<Nfz ; ++n ){
      face_riemann_z( face_pointer(theFaces_z,n) , dt );
    }
  }
  //Source Terms
  cell_add_src( theCells ,theGrid, theGravMasses , dt );
  //forceGravMasss( theCells , theGravMasses );

  //Bookkeeping
  cell_update_phi( theCells ,theGrid, RK , dt );
  cell_update_dphi( theCells ,theGrid);
  // update_RK_gravMasses( theGravMasses , RK , dt );
  cell_calc_prim( theCells ,theGrid);

  //Boundary Data
  //cell_boundary_outflow_r( theCells , theFaces_r ,theGrid, nri );
  //if( N_z_global > 1 ) cell_boundary_z( theCells , theFaces_z ,theGrid, nzk );
  cell_boundary_fixed_r( theCells, theGrid );

  free(nri);
  free(nzk);
  face_destroy(theFaces_r);
  face_destroy(theFaces_z);

  //inter-processor syncs
  cell_syncproc_r(theCells,theGrid);
  cell_syncproc_z(theCells,theGrid);
  
  cell_prim2cons( theCells,theGrid );


}

void onestep_Psi( struct Cell *** theCells , struct Grid * theGrid, double dt ){
  int N_r_withghost = grid_N_r(theGrid)+grid_Nghost_rmin(theGrid)+grid_Nghost_rmax(theGrid);
  int N_z_withghost = grid_N_z(theGrid)+grid_Nghost_zmin(theGrid)+grid_Nghost_zmax(theGrid);

  int i,j,k;
  for( k=0 ; k<N_z_withghost ; ++k ){
    for( i=0 ; i<N_r_withghost ; ++i ){
      for( j=0 ; j<grid_N_p(theGrid,i) ; ++j ){
        cell_mult_psi(cell_pointer(theCells,i,j,k),exp(-dt*DIVB_CH/DIVB_L));
      }
    }
  }
  //Boundary Data: WE MAY NEED TO APPLY BCs on PSI NOW.

  //Bookkeeping
  cell_calc_prim( theCells,theGrid );
  cell_prim2cons( theCells,theGrid );

  //inter-processor syncs
  cell_syncproc_r(theCells,theGrid);
  cell_syncproc_z(theCells,theGrid);

}


