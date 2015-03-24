#define TIMESTEP_PRIVATE_DEFS
#define RIEMANN_PRIVATE_DEFS
#define CELL_PRIVATE_DEFS
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

  int i,j,k;
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

  //TODO: Remove these dumb print statements
  if(PRINTTOOMUCH)
  {
    cell_print_all(theCells, theSim);
  }
  
  // this is part of the runge-kutta method
  cell_adjust_RK_cons( theCells, theSim, theTimeStep->RK);
  
  //TODO: Remove these dumb print statements
  
  if(PRINTTOOMUCH)
  {
    cell_print_all(theCells, theSim);
  }
  

  //Piecewise-Linear Reconstructions.
  //  must be done before AddFlux to find gradients for viscosity
  cell_plm_p(theCells,theSim);
  cell_plm_rz(theCells,theSim,theFaces_r,theTimeStep,theMPIsetup,R_DIR);
  if( sim_N_global(theSim,Z_DIR) != 1 )
    cell_plm_rz(theCells,theSim,theFaces_z,theTimeStep,theMPIsetup,Z_DIR );

  if(PRINTTOOMUCH && 0)
  {
    FILE *gradfile = fopen("grad.out", "w");
    for( k=0 ; k<sim_N(theSim,Z_DIR) ; ++k ){
      double zp = sim_FacePos(theSim,k,Z_DIR);
      double zm = sim_FacePos(theSim,k-1,Z_DIR);
      double z = 0.5*(zp+zm);
      for( i=0 ; i<sim_N(theSim,R_DIR) ; ++i ){
        double rp = sim_FacePos(theSim,i,R_DIR);
        double rm = sim_FacePos(theSim,i-1,R_DIR);
        double r = 0.5*(rp+rm);
        for( j=0 ; j<sim_N_p(theSim,i) ; ++j ){   
          struct Cell * c = &(theCells[k][i][j]);
          double p = cell_tiph(c) - 0.5*cell_dphi(c);
          fprintf(gradfile, "%d, %d, %d, %.12g, %.12g, %.12g, %.12g, %.12g, %.12g, %.12g, %.12g, %.12g, %.12g, %.12g, %.12g, %.12g, %.12g, %.12g, %.12g, %.12g, %.12g, %.12g, %.12g, %.12g, %.12g, %.12g\n",i,j,k,r,p,z,c->prim[RHO],c->prim[URR],c->prim[UPP],c->prim[UZZ],c->prim[PPP],c->gradr[RHO],c->gradr[URR],c->gradr[UPP],c->gradr[UZZ],c->gradr[PPP],c->gradp[RHO],c->gradp[URR],c->gradp[UPP],c->gradp[UZZ],c->gradp[PPP],c->gradz[RHO],c->gradz[URR],c->gradz[UPP],c->gradz[UZZ],c->gradz[PPP]);
          }
      }
    }  
    fclose(gradfile);
    gradfile = fopen("grad_face.out","w");
    fclose(gradfile);
  }

  //Phi Flux
  for( k=0 ; k<sim_N(theSim,Z_DIR) ; ++k ){
    for( i=0 ; i<sim_N(theSim,R_DIR) ; ++i ){
      for( j=0 ; j<sim_N_p(theSim,i) ; ++j ){
        struct Riemann * theRiemann = riemann_create(theSim); // struct to contain everything we need to solve Riemann problem 
        riemann_setup_p(theRiemann,theCells,theSim,i,j,k,PDIRECTION); // set various quantities in theRiemann
        riemann_AddFlux(theRiemann,theSim,dt); // solve Riemann problem and update RHS
       
        if(PRINTTOOMUCH)
        {
            double *pL = theRiemann->primL;
            double *pR = theRiemann->primR;
            double *FL = theRiemann->FL;
            double *FR = theRiemann->FR;
            double *us = theRiemann->Ustar;
            double *fs = theRiemann->F;
            printf("Face: r=%.12g, phi=%.12g\n",theRiemann->pos[R_DIR],0.5*(cell_tiph(theRiemann->cL)-cell_dphi(theRiemann->cL)+cell_tiph(theRiemann->cR)-cell_dphi(theRiemann->cR)));
            printf("Flux Phi L (rho=%.12g, Pp=%.12g, vr=%.12g, vp=%.12g, vz=%.12g): DDD:%.12g, SRR:%.12g, LLL:%.12g, SZZ:%.12g, TAU:%.12g\n", pL[RHO], pL[PPP], pL[URR], pL[UPP], pL[UZZ],FL[DDD],FL[SRR],FL[LLL],FL[SZZ],FL[TAU]);
            printf("Flux Phi R (rho=%.12g, Pp=%.12g, vr=%.12g, vp=%.12g, vz=%.12g): DDD:%.12g, SRR:%.12g, LLL:%.12g, SZZ:%.12g, TAU:%.12g\n", pR[RHO], pR[PPP], pR[URR], pR[UPP],pR[UZZ],FR[DDD],FR[SRR],FR[LLL],FR[SZZ],FR[TAU]);
            printf("Flux Phi * (rhostar=%.12g, Sr=%.12g, Sp=%.12g, Sz=%.12g, tau=%.12g): DDD:%.12g, SRR:%.12g, LLL:%.12g, SZZ:%.12g, TAU:%.12g\n", us[DDD], us[SRR], us[LLL], us[SZZ], us[TAU],fs[DDD],fs[SRR],fs[LLL],fs[SZZ],fs[TAU]);
        if(fs[SZZ] != 0.0)
            printf("   SZZ flux!\n");
        if(us[SZZ] != 0.0)
            printf("   SZZ on face!\n");
        }
        riemann_destroy(theRiemann); // clean up

      }
    }
  }

  //R Flux
  int n;
  for( n=0 ; n<timestep_n(theTimeStep,sim_N(theSim,R_DIR)-1,R_DIR) ; ++n ){
    struct Riemann * theRiemann = riemann_create(theSim); //struct to contain everything we need to solve Riemann problem
    riemann_setup_rz(theRiemann,theFaces_r,theSim,n,RDIRECTION);  //set various quantities in theRiemann
    riemann_AddFlux(theRiemann,theSim,dt); // solve Riemann problem and update RHS
    
    if(PRINTTOOMUCH)
    {
        double *pL = theRiemann->primL;
        double *pR = theRiemann->primR;
        double *FL = theRiemann->FL;
        double *FR = theRiemann->FR;
        double *us = theRiemann->Ustar;
        double *fs = theRiemann->F;
        double Sl = theRiemann->Sl;
        double Sr = theRiemann->Sr;
        printf("Face: r=%.12g, phi=%.12g, dA=%.12g, cm=%.12g\n",theRiemann->pos[R_DIR],0.5*(cell_tiph(theRiemann->cL)-0.5*cell_dphi(theRiemann->cL)+cell_tiph(theRiemann->cR)-0.5*cell_dphi(theRiemann->cR)), face_dA(theFaces_r,n), face_cm(theFaces_r, n));
        printf("Flux Rad L (Sl=%.12g) (rho=%.12g, Pp=%.12g, vr=%.12g, vp=%.12g, vz=%.12g): DDD:%.12g, SRR:%.12g, LLL:%.12g, SZZ:%.12g, TAU:%.12g\n", Sl, pL[RHO], pL[PPP], pL[URR], pL[UPP],pL[UZZ],FL[DDD],FL[SRR],FL[LLL],FL[SZZ],FL[TAU]);
        printf("Flux Rad R (Sr=%.12g) (rho=%.12g, Pp=%.12g, vr=%.12g, vp=%.12g, vz=%.12g): DDD:%.12g, SRR:%.12g, LLL:%.12g, SZZ:%.12g, TAU:%.12g\n", Sr, pR[RHO], pR[PPP], pR[URR], pR[UPP],pR[UZZ],FR[DDD],FR[SRR],FR[LLL],FR[SZZ],FR[TAU]);
        printf("Flux Rad * (rhostar=%.12g, Sr=%.12g, Sp=%.12g, Sz=%.12g, tau=%.12g): DDD:%.12g, SRR:%.12g, LLL:%.12g, SZZ:%.12g, TAU:%.12g\n", us[DDD], us[SRR], us[LLL], us[SZZ], us[TAU],fs[DDD],fs[SRR],fs[LLL],fs[SZZ],fs[TAU]);
    }
    
    riemann_destroy(theRiemann); // clean up
  }
  //Z Flux
  if( sim_N_global(theSim,Z_DIR) != 1 ){
      if(PRINTTOOMUCH)
          printf("WHOA, Z FLUX HAPPENING!\n");
    for( n=0 ; n<timestep_n(theTimeStep,sim_N(theSim,Z_DIR)-1,Z_DIR); ++n ){
      struct Riemann * theRiemann = riemann_create(theSim); // struct to contain everything we need to solve Riemann problem
      riemann_setup_rz(theRiemann,theFaces_z,theSim,n,ZDIRECTION); // set various quantities in theRiemann
      riemann_AddFlux(theRiemann,theSim,dt); // solve Riemann problem and update RHS 
      riemann_destroy(theRiemann); // clean up
    }
  }
  
  //TODO: Remove these dumb print statements
  if(PRINTTOOMUCH)
  {
    cell_print_all(theCells, theSim);
  }
  
  //Source Terms
  cell_add_src( theCells ,theSim, theGravMasses , dt ); // add source terms

    
  //TODO: Remove these dumb print statements
  if(PRINTTOOMUCH)
  {
    cell_print_all(theCells, theSim);
  }
  

  //Bookkeeping
  cell_update_phi( theCells ,theSim, theTimeStep->RK , dt ); // allow the cells to move in phi direction
  cell_update_dphi( theCells ,theSim); // all the cells to change size
  gravMass_update_RK( theGravMasses ,theSim, theTimeStep->RK ); // allow the GravMasses to move
  cell_calc_prim( theCells ,theSim); // calculate primitives

  //inter-processor syncs
  cell_syncproc_r(theCells,theSim,theMPIsetup);
  cell_syncproc_z(theCells,theSim,theMPIsetup);

  //Boundary Data
  //R - Inner
  if (sim_BoundTypeRIn(theSim)==BOUND_OUTFLOW)
    cell_boundary_outflow_r_inner( theCells, theFaces_r, theSim, theMPIsetup, 
                                    theTimeStep);
  else if (sim_BoundTypeRIn(theSim)==BOUND_FIXED)
    cell_boundary_fixed_r_inner(theCells, theSim, theMPIsetup, 
                                    (*cell_single_init_ptr(theSim)));
  else if (sim_BoundTypeRIn(theSim)==BOUND_SS)
    cell_boundary_ssprofile_r_inner(theCells, theSim, theMPIsetup);    
  else if (sim_BoundTypeRIn(theSim)==BOUND_LINEAR)
    cell_boundary_linear_r_inner( theCells, theFaces_r, theSim, theMPIsetup, 
                                    theTimeStep);
  
  //R - Outer
  if (sim_BoundTypeROut(theSim)==BOUND_OUTFLOW)
    cell_boundary_outflow_r_outer( theCells, theFaces_r, theSim, theMPIsetup, 
                                    theTimeStep );
  else if (sim_BoundTypeROut(theSim)==BOUND_FIXED)
    cell_boundary_fixed_r_outer( theCells, theSim, theMPIsetup,
                                (*cell_single_init_ptr(theSim)));
  else if (sim_BoundTypeROut(theSim)==BOUND_SS)
    cell_boundary_ssprofile_r_outer(theCells, theSim, theMPIsetup);
  else if (sim_BoundTypeROut(theSim)==BOUND_LINEAR)
    cell_boundary_linear_r_outer( theCells, theFaces_r, theSim, theMPIsetup, 
                                    theTimeStep);
  
  if (sim_N_global(theSim,Z_DIR)>1){
    //Z - Bottom
    if (sim_BoundTypeZBot(theSim)==BOUND_OUTFLOW){
      cell_boundary_outflow_z_bot( theCells , theFaces_z ,theSim,theMPIsetup, theTimeStep );
    }else if (sim_BoundTypeZBot(theSim)==BOUND_FIXED){
      cell_boundary_fixed_z_bot(theCells,theSim,theMPIsetup,(*cell_single_init_ptr(theSim)));    
    } else if (sim_BoundTypeZBot(theSim)==BOUND_PERIODIC){
      //do nothing, this is already handled by the syncing routine
    }
    //Z - Top
    if (sim_BoundTypeZTop(theSim)==BOUND_OUTFLOW){
      cell_boundary_outflow_z_top( theCells , theFaces_z ,theSim,theMPIsetup, theTimeStep );
    }else if (sim_BoundTypeZTop(theSim)==BOUND_FIXED){
      cell_boundary_fixed_z_top(theCells,theSim,theMPIsetup,(*cell_single_init_ptr(theSim)));    
    } else if (sim_BoundTypeZTop(theSim)==BOUND_PERIODIC){
      //do nothing, this is already handled by the syncing routine
    }
  }

  // Add External Sources
  if(sim_BoundTypeSource(theSim) == BOUND_NOZZLE)
      cell_boundary_nozzle(theCells, theSim, theMPIsetup, theTimeStep);

  //clean up
  face_destroy(theFaces_r);
  if (sim_N_global(theSim,Z_DIR)>1){
    face_destroy(theFaces_z);
  }

  // if DAMP_TIME is set to a positive number, apply damping near boundary
  if (sim_DAMP_TIME(theSim)>0.0) cell_bc_damp( theCells , theSim, dt,(*cell_single_init_ptr(theSim)) );

  //re-calculate conserved quantities. 
  //things may have changed due to syncing and/or caps/floors applied in primitive solver.
  cell_calc_cons( theCells,theSim );
}


