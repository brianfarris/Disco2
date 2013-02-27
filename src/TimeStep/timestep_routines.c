#define TIMESTEP_PRIVATE_DEFS
#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/TimeStep.h"
#include "../Headers/Grid.h"
#include "../Headers/Cell.h"
#include "../Headers/GravMass.h"
#include "../Headers/Face.h"
#include "../Headers/TimeStep.h"
#include "../Headers/MPIsetup.h"
#include "../Headers/header.h"

void timestep_set_dt(struct TimeStep * theTimeStep, struct Cell *** theCells, struct Grid * theGrid){
  theTimeStep->dt = cell_mindt(theCells,theGrid);
  if( theTimeStep->t+theTimeStep->dt > theTimeStep->T_MAX ) {
    theTimeStep->dt = theTimeStep->T_MAX-theTimeStep->t;
  }
  printf("t: %e, dt: %e\n",theTimeStep->t,theTimeStep->dt);
}
void timestep_update_t(struct TimeStep * theTimeStep){
  theTimeStep->t += theTimeStep->dt;
}
void timestep_set_RK(struct TimeStep * theTimeStep,double RK){
  theTimeStep->RK = RK;
}
double timestep_get_t(struct TimeStep * theTimeStep){
  return(theTimeStep->t);
}
double timestep_get_T_MAX(struct TimeStep * theTimeStep){
  return(theTimeStep->T_MAX);
}
double timestep_NUM_CHECKPOINTS(struct TimeStep * theTimeStep){
  return(theTimeStep->NUM_CHECKPOINTS);
}



void timestep_substep(struct TimeStep * theTimeStep, struct Cell *** theCells,struct Grid * theGrid,struct GravMass * theGravMasses,struct MPIsetup * theMPIsetup,double timestep_fac){

  double dt = timestep_fac*theTimeStep->dt;
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
  if( grid_MOVE_CELLS(theGrid) == C_WCELL ) cell_set_wcell( theCells ,theGrid);
  if( grid_MOVE_CELLS(theGrid) == C_RIGID ) cell_set_wrigid( theCells ,theGrid);
  cell_adjust_RK_cons( theCells, theGrid, theTimeStep->RK);
  cell_clear_divB(theCells,theGrid);
  cell_clear_GradPsi(theCells,theGrid);
  //Phi Flux
  cell_plm_p( theCells,theGrid );
  cell_flux_p( theCells ,theGrid, dt );
  //R Flux
  cell_plm_rz( theCells ,theGrid, theFaces_r , Nfr , 0 );
  int n;

  for( n=0 ; n<Nfr ; ++n ){
    face_riemann_r( face_pointer(theFaces_r,n) , theGrid, dt );
  }
  //Z Flux
  if( grid_N_z_global(theGrid) != 1 ){
    cell_plm_rz( theCells ,theGrid, theFaces_z , Nfz , 1 );
    for( n=0 ; n<Nfz ; ++n ){
      face_riemann_z( face_pointer(theFaces_z,n) ,theGrid, dt );
    }
  }
  //Source Terms
  cell_add_src( theCells ,theGrid, theGravMasses , dt );
  //forceGravMasss( theCells , theGravMasses );

  //Bookkeeping
  cell_update_phi( theCells ,theGrid, theTimeStep->RK , dt );
  cell_update_dphi( theCells ,theGrid);
  // update_RK_gravMasses( theGravMasses , RK , dt );
  cell_calc_prim( theCells ,theGrid);

  //Boundary Data
  //cell_boundary_outflow_r( theCells , theFaces_r ,theGrid, nri );
  //if( N_z_global > 1 ) cell_boundary_z( theCells , theFaces_z ,theGrid, nzk );
  cell_boundary_fixed_r( theCells, theGrid,theMPIsetup );

  free(nri);
  free(nzk);
  face_destroy(theFaces_r);
  face_destroy(theFaces_z);

  //inter-processor syncs
  cell_syncproc_r(theCells,theGrid,theMPIsetup);
  cell_syncproc_z(theCells,theGrid,theMPIsetup);

  cell_calc_cons( theCells,theGrid );


}

void timestep_update_Psi( struct TimeStep * theTimeStep, struct Cell *** theCells , struct Grid * theGrid,struct MPIsetup * theMPIsetup){
  int N_r_withghost = grid_N_r(theGrid)+grid_Nghost_rmin(theGrid)+grid_Nghost_rmax(theGrid);
  int N_z_withghost = grid_N_z(theGrid)+grid_Nghost_zmin(theGrid)+grid_Nghost_zmax(theGrid);
  double DIVB_CH = grid_DIVB_CH(theGrid);
  double DIVB_L = grid_DIVB_L(theGrid);

  int i,j,k;
  for( k=0 ; k<N_z_withghost ; ++k ){
    for( i=0 ; i<N_r_withghost ; ++i ){
      for( j=0 ; j<grid_N_p(theGrid,i) ; ++j ){
        cell_mult_psi(cell_pointer(theCells,i,j,k),exp(-theTimeStep->dt*DIVB_CH/DIVB_L));
      }
    }
  }
  //Boundary Data: WE MAY NEED TO APPLY BCs on PSI NOW.

  //Bookkeeping
  cell_calc_prim( theCells,theGrid );
  cell_calc_cons( theCells,theGrid );

  //inter-processor syncs
  cell_syncproc_r(theCells,theGrid,theMPIsetup);
  cell_syncproc_z(theCells,theGrid,theMPIsetup);

}


