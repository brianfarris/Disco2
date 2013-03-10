#define TIMESTEP_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/TimeStep.h"
#include "../Headers/Riemann.h"
#include "../Headers/Sim.h"
#include "../Headers/Cell.h"
#include "../Headers/GravMass.h"
#include "../Headers/Face.h"
#include "../Headers/TimeStep.h"
#include "../Headers/MPIsetup.h"
#include "../Headers/header.h"

void timestep_set_dt(struct TimeStep * theTimeStep, struct Cell *** theCells, struct Sim * theSim){
  theTimeStep->dt = cell_mindt(theCells,theSim);
  if( theTimeStep->t+theTimeStep->dt > sim_get_T_MAX(theSim) ) {
    theTimeStep->dt = sim_get_T_MAX(theSim)-theTimeStep->t;
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
double timestep_dt(struct TimeStep * theTimeStep){
  return(theTimeStep->dt);
}

void timestep_substep(struct TimeStep * theTimeStep, struct Cell *** theCells,struct Sim * theSim,struct GravMass * theGravMasses,struct MPIsetup * theMPIsetup,double timestep_fac){

  double dt = timestep_fac*theTimeStep->dt;
  struct Face *theFaces_r = face_create(theCells,theSim,theTimeStep,0); // r-direction
  struct Face *theFaces_z = face_create(theCells,theSim,theTimeStep,1); // z-direction
  cell_clean_pi(theCells,theSim);
  gravMass_clean_pi(theGravMasses,theSim);
  cell_clear_w(theCells,theSim);
  if( sim_MOVE_CELLS(theSim) == C_WCELL ) cell_set_wcell( theCells ,theSim);
  if( sim_MOVE_CELLS(theSim) == C_RIGID ) cell_set_wrigid( theCells ,theSim);
  cell_adjust_RK_cons( theCells, theSim, theTimeStep->RK);
  cell_clear_divB(theCells,theSim);
  cell_clear_GradPsi(theCells,theSim);
  //Phi Flux
  cell_plm_p( theCells,theSim );
  int i,j,k;
  for( k=0 ; k<sim_N(theSim,Z_DIR) ; ++k ){
    for( i=0 ; i<sim_N(theSim,R_DIR) ; ++i ){
      for( j=0 ; j<sim_N_p(theSim,i)-1 ; ++j ){
        struct Riemann * theRiemann = riemann_create(theSim);
        riemann_setup_p(theRiemann,theCells,theSim,i,j,j+1,k);
        riemann_hllc( theRiemann,theSim,dt,1);
        riemann_destroy(theRiemann);
      }
      struct Riemann * theRiemann = riemann_create(theSim);
      riemann_setup_p(theRiemann,theCells,theSim,i,sim_N_p(theSim,i)-1,0,k);
      riemann_hllc(theRiemann, theSim,dt,1 );
      riemann_destroy(theRiemann);
    }
  }
  //R Flux
  cell_plm_rz( theCells ,theSim, theFaces_r , timestep_Nfr(theTimeStep) , 0 );
  int n;
  for( n=0 ; n<timestep_Nfr(theTimeStep) ; ++n ){
    struct Riemann * theRiemann = riemann_create(theSim);
    riemann_setup_rz(theRiemann,theFaces_r,theSim,n);
    riemann_hllc(theRiemann,theSim,dt,0);
    riemann_destroy(theRiemann);
  }

  //Z Flux
  if( sim_N_global(theSim,Z_DIR) != 1 ){
    cell_plm_rz( theCells ,theSim, theFaces_z , timestep_Nfz(theTimeStep) , 1 );
    for( n=0 ; n<timestep_Nfz(theTimeStep) ; ++n ){
      struct Riemann * theRiemann = riemann_create(theSim);
      riemann_setup_rz(theRiemann,theFaces_z,theSim,n);
      riemann_hllc(theRiemann,theSim,dt,2); 
      riemann_destroy(theRiemann);
    }
  }

  //Source Terms
  cell_add_src( theCells ,theSim, theGravMasses , dt );
  if (sim_EXPLICIT_VISCOSITY(theSim)>0.0){
    cell_add_visc_src( theCells ,theSim,dt );
  }

  //Bookkeeping
  cell_update_phi( theCells ,theSim, theTimeStep->RK , dt );
  cell_update_dphi( theCells ,theSim);
  gravMass_update_RK( theGravMasses ,theSim, theTimeStep->RK );
  cell_calc_prim( theCells ,theSim);

  //Boundary Data
  //if( N_z_global > 1 ) cell_boundary_z( theCells , theFaces_z ,theSim, nzk );
  if (sim_BoundTypeR(theSim)==BOUND_OUTFLOW){
    cell_boundary_outflow_r( theCells , theFaces_r ,theSim,theMPIsetup, theTimeStep->nri );
  }else if (sim_BoundTypeR(theSim)==BOUND_FIXED){
    cell_boundary_fixed_r(theCells,theSim,theMPIsetup,(*cell_single_init_ptr(theSim)));    
  }
  if (sim_BoundTypeZ(theSim)==BOUND_OUTFLOW){
    cell_boundary_outflow_z( theCells , theFaces_r ,theSim,theMPIsetup, theTimeStep->nzk );
  }else if (sim_BoundTypeZ(theSim)==BOUND_FIXED){
    cell_boundary_fixed_z(theCells,theSim,theMPIsetup,(*cell_single_init_ptr(theSim)));    
  } else if (sim_BoundTypeZ(theSim)==BOUND_PERIODIC){
    //do nothing, this is already handled by the syncing routine
  }
   

  face_destroy(theFaces_r);
  if (sim_N_global(theSim,Z_DIR)>1){
    face_destroy(theFaces_z);
  }

  //inter-processor syncs
  cell_syncproc_r(theCells,theSim,theMPIsetup);
  if (sim_N_global(theSim,Z_DIR)>1){
    cell_syncproc_z(theCells,theSim,theMPIsetup);
  }
  
  cell_calc_cons( theCells,theSim );

}

void timestep_update_Psi( struct TimeStep * theTimeStep, struct Cell *** theCells , struct Sim * theSim,struct MPIsetup * theMPIsetup){
  double DIVB_CH = sim_DIVB_CH(theSim);
  double DIVB_L = sim_DIVB_L(theSim);

  int i,j,k;
  for( k=0 ; k<sim_N(theSim,Z_DIR) ; ++k ){
    for( i=0 ; i<sim_N(theSim,R_DIR) ; ++i ){
      for( j=0 ; j<sim_N_p(theSim,i) ; ++j ){
        cell_mult_psi(cell_single(theCells,i,j,k),exp(-theTimeStep->dt*DIVB_CH/DIVB_L));
      }
    }
  }
  //Boundary Data: WE MAY NEED TO APPLY BCs on PSI NOW.

  //Bookkeeping
  cell_calc_prim( theCells,theSim );
  cell_calc_cons( theCells,theSim );

  //inter-processor syncs
  cell_syncproc_r(theCells,theSim,theMPIsetup);
  if (sim_N_global(theSim,Z_DIR)>1){
    cell_syncproc_z(theCells,theSim,theMPIsetup);
  }
}

int * timestep_nri(struct TimeStep * theTimeStep){
  return(theTimeStep->nri);
}
int * timestep_nzk(struct TimeStep * theTimeStep){
  return(theTimeStep->nzk);
}

void timestep_set_Nfr(struct TimeStep * theTimeStep,struct Sim * theSim){
  theTimeStep->Nfr = theTimeStep->nri[sim_N(theSim,R_DIR)-1];
}
void timestep_set_Nfz(struct TimeStep * theTimeStep,struct Sim * theSim){
  theTimeStep->Nfz = theTimeStep->nzk[sim_N(theSim,Z_DIR)-1];
}

int timestep_Nfr(struct TimeStep * theTimeStep){
  return(theTimeStep->Nfr);
}
int timestep_Nfz(struct TimeStep * theTimeStep){
  return(theTimeStep->Nfz);
}


