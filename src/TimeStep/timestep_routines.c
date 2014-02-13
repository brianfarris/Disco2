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

void timestep_substep(struct TimeStep * theTimeStep, struct Cell *** theCells,
    struct Sim * theSim,struct GravMass * theGravMasses,
    struct MPIsetup * theMPIsetup,double timestep_fac){

  double dt = timestep_fac*theTimeStep->dt;

  // figure out what all the faces need to be and set them up
  struct Face *theFaces_r = face_create(theCells,theSim,theTimeStep,0); // r-direction
  struct Face *theFaces_z = face_create(theCells,theSim,theTimeStep,1); // z-direction

  // make sure phi is between 0 and 2pi for cells and GravMasses. Not sure why we care.
  cell_clean_pi(theCells,theSim);
  gravMass_clean_pi(theGravMasses,theSim);

  // reset wiph
  cell_clear_w(theCells,theSim);
  cell_set_w( theCells ,theSim);

  // this is part of the runge-kutta method
  cell_adjust_RK_cons( theCells, theSim, theTimeStep->RK);

  // divB and GradPsi are needed for Powell source terms in MHD eqns. Here we reset them to 0.
  cell_clear_divB(theCells,theSim);
  cell_clear_GradPsi(theCells,theSim);

  double r0 = sim_FacePos(theSim,0,R_DIR);
  double r1 = sim_FacePos(theSim,1,R_DIR);
  double r2 = sim_FacePos(theSim,2,R_DIR);

  int i,j,k,q;
  //Phi Flux
  cell_plm_p(theCells,theSim);//piecewise-linear reconstruction
  for( k=0 ; k<sim_N(theSim,Z_DIR) ; ++k ){
    for( i=0 ; i<sim_N(theSim,R_DIR) ; ++i ){
      for( j=0 ; j<sim_N_p(theSim,i) ; ++j ){
        struct Riemann * theRiemann = riemann_create(theSim); // struct to contain everything we need to solve Riemann problem 
        riemann_setup_p(theRiemann,theCells,theSim,i,j,k,PDIRECTION); // set various quantities in theRiemann
        riemann_AddFlux(theRiemann,theSim,dt); // solve Riemann problem and update RHS
        riemann_destroy(theRiemann); // clean up
      }
    }
  }

  //R Flux
  cell_plm_rz(theCells,theSim,theFaces_r,theTimeStep,theMPIsetup,R_DIR); //piecewise-linear reconstruction

  int n;
  for( n=0 ; n<timestep_n(theTimeStep,sim_N(theSim,R_DIR)-1,R_DIR) ; ++n ){
    struct Riemann * theRiemann = riemann_create(theSim); //struct to contain everything we need to solve Riemann problem
    riemann_setup_rz(theRiemann,theFaces_r,theSim,n,RDIRECTION);  //set various quantities in theRiemann
    riemann_AddFlux(theRiemann,theSim,dt); // solve Riemann problem and update RHS
    riemann_destroy(theRiemann); // clean up
  }

  //Z Flux
  if( sim_N_global(theSim,Z_DIR) != 1 ){
    cell_plm_rz(theCells,theSim,theFaces_z,theTimeStep,theMPIsetup,Z_DIR );//piecewise-linear reconstruction
    for( n=0 ; n<timestep_n(theTimeStep,sim_N(theSim,Z_DIR)-1,Z_DIR); ++n ){
      struct Riemann * theRiemann = riemann_create(theSim); // struct to contain everything we need to solve Riemann problem
      riemann_setup_rz(theRiemann,theFaces_z,theSim,n,ZDIRECTION); // set various quantities in theRiemann
      riemann_AddFlux(theRiemann,theSim,dt); // solve Riemann problem and update RHS 
      riemann_destroy(theRiemann); // clean up
    }
  }

  //Source Terms
  cell_add_src( theCells ,theSim, theGravMasses , dt ); // add source terms
  if (sim_EXPLICIT_VISCOSITY(theSim)>0.0){
    if (VISC_OLD==1){
      cell_add_visc_src_old( theCells ,theSim,dt ); // add viscous source terms
    } else{
      cell_add_visc_src( theCells ,theSim,dt ); // add viscous source terms
    }
  }

  //cell_add_split_fictitious(theCells,theSim,dt);

  //Bookkeeping
  cell_update_phi( theCells ,theSim, theTimeStep->RK , dt ); // allow the cells to move in phi direction
  cell_update_dphi( theCells ,theSim); // all the cells to change size
  gravMass_update_RK( theGravMasses ,theSim, theTimeStep->RK ); // allow the GravMasses to move
  cell_calc_prim( theCells ,theSim); // calculate primitives

  //Set temperature to keep HoR constant
  if (sim_SET_T(theSim)==1){
    cell_setT(theCells,theSim,theGravMasses); 
  }

  //Boundary Data
  if (sim_BoundTypeR(theSim)==BOUND_OUTFLOW){
    cell_boundary_outflow_r( theCells , theFaces_r ,theSim,theMPIsetup, theTimeStep );
  }else if (sim_BoundTypeR(theSim)==BOUND_FIXED){
    cell_boundary_fixed_r(theCells,theSim,theMPIsetup,(*cell_single_init_ptr(theSim)));    
  }
  if (sim_N_global(theSim,Z_DIR)>1){
    if (sim_BoundTypeZ(theSim)==BOUND_OUTFLOW){
      cell_boundary_outflow_z( theCells , theFaces_z ,theSim,theMPIsetup, theTimeStep );
    }else if (sim_BoundTypeZ(theSim)==BOUND_FIXED){
      cell_boundary_fixed_z(theCells,theSim,theMPIsetup,(*cell_single_init_ptr(theSim)));    
    } else if (sim_BoundTypeZ(theSim)==BOUND_PERIODIC){
      //do nothing, this is already handled by the syncing routine
    }
  } 
  
  //clean up
  face_destroy(theFaces_r);
  if (sim_N_global(theSim,Z_DIR)>1){
    face_destroy(theFaces_z);
  }

  // if DAMP_TIME is set to a positive number, apply damping near boundary
  if (sim_DAMP_TIME(theSim)>0.0) cell_bc_damp( theCells , theSim, dt,(*cell_single_init_ptr(theSim)) );

  //inter-processor syncs
  cell_syncproc_r(theCells,theSim,theMPIsetup);
  cell_syncproc_z(theCells,theSim,theMPIsetup);

  //re-calculate conserved quantities. 
  //things may have changed due to syncing and/or caps/floors applied in primitive solver.
  
  cell_calc_cons( theCells,theSim );

}

void timestep_update_Psi( struct TimeStep * theTimeStep, struct Cell *** theCells , struct Sim * theSim,struct MPIsetup * theMPIsetup){
  double DIVB_CH = sim_DIVB_CH(theSim);
  double DIVB_L = sim_DIVB_L(theSim);

  int i,j,k;
  for( k=0 ; k<sim_N(theSim,Z_DIR) ; ++k ){
    for( i=0 ; i<sim_N(theSim,R_DIR) ; ++i ){
      for( j=0 ; j<sim_N_p(theSim,i) ; ++j ){
        // we are handling Psi equation in an operator splitting manner. 
        // For this step the solution is analytic. See Dedner paper for discussion.
        cell_mult_psi(cell_single(theCells,i,j,k),exp(-theTimeStep->dt*DIVB_CH/DIVB_L)); 
      }
    }
  }

  //SHOULD WE APPLY BCs on PSI NOW?

  //Bookkeeping. Is this really necessary?
  cell_calc_prim( theCells,theSim );
  cell_calc_cons( theCells,theSim );

  //inter-processor syncs
  cell_syncproc_r(theCells,theSim,theMPIsetup);
  if (sim_N_global(theSim,Z_DIR)>1){
    cell_syncproc_z(theCells,theSim,theMPIsetup);
  }
}


