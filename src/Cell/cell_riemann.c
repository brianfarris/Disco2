#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Grid.h"
#include "../Headers/Flux.h"
#include "../Headers/GravMass.h"
#include "../Headers/header.h"

void vel( double * , double * , double * , double * , double * , double * , double , double *,double ,double );
void getUstar( double * , double * , double , double , double , double * , double * ,double);
void flux( double * , double * , double , double * ,double ,double );

void cell_riemann_p( struct Cell * cL , struct Cell * cR, struct Grid * theGrid, double dA , double dt , double r ){
  int NUM_Q = grid_NUM_Q(theGrid);
  double GAMMALAW = grid_GAMMALAW(theGrid);
  double DIVB_CH = grid_DIVB_CH(theGrid);

  struct Flux * theFlux = flux_create(theGrid);
/*
  double * primL = malloc(NUM_Q*sizeof(double));
  double * primR = malloc(NUM_Q*sizeof(double));

  int q;
  for( q=0 ; q<NUM_Q ; ++q ){
    primL[q] = cell_prims(cL)[q] 
      + cell_gradp(cL)[q]*cell_dphi(cL);
    primR[q] = cell_prims(cR)[q] 
      + cell_gradp(cR)[q]*cell_dphi(cR);
  }
  */

  int q;
  for (q=0;q<NUM_Q;++q){
    flux_set_primL(theFlux,q,cell_prims(cL)[q] + cell_gradp(cL)[q]*cell_dphi(cL));
    flux_set_primR(theFlux,q,cell_prims(cR)[q] + cell_gradp(cR)[q]*cell_dphi(cR));
  }

  double Sl,Sr,Ss;
  double n[3];
  n[0] = 0.0;   n[1] = 1.0;   n[2] = 0.0;

  double Bpack[6];
  //vel( primL , primR , &Sl , &Sr , &Ss , n , r , Bpack ,GAMMALAW,DIVB_CH);
  flux_set_vel(theFlux,n,r,Bpack,GAMMALAW,DIVB_CH);

  //double * Fk = malloc(NUM_Q*sizeof(double));
  double * Uk = malloc(NUM_Q*sizeof(double));

  //double * Flux = malloc(NUM_Q*sizeof(double));


  double Bp_face = 0.5*(flux_primL(theFlux)[BPP]+flux_primR(theFlux)[BPP]);
  double Psi_face = 0.5*(flux_primL(theFlux)[PSI]+flux_primR(theFlux)[PSI]);

  if( grid_MOVE_CELLS(theGrid) == C_WRIEMANN ) cell_add_wiph(cL,Ss);
  double w = cell_wiph(cL);
  if( w < Sl ){
    flux( primL , Fk , r , n ,GAMMALAW,DIVB_CH);
    cell_prim2cons( primL , Uk , r , 1.0,GAMMALAW );
    flux_set_Uk(theFlux,Uk);
    for( q=0 ; q<NUM_Q ; ++q ){
      Flux[q] = Fk[q] - w*Uk[q];
    }
  }else if( w > Sr ){
    flux( primR , Fk , r , n,GAMMALAW,DIVB_CH );
    cell_prim2cons( primR , Uk , r , 1.0 ,GAMMALAW);
    flux_set_Uk(theFlux,Uk);
    for( q=0 ; q<NUM_Q ; ++q ){
      Flux[q] = Fk[q] - w*Uk[q];
    }
  }else{
    if( w < Ss ){
      cell_prim2cons( primL , Uk , r , 1.0 ,GAMMALAW);
    }else{
      cell_prim2cons( primR , Uk , r , 1.0 ,GAMMALAW);        
    }
    flux_set_Uk(theFlux,Uk);
    flux_set_Ustar(theFlux,n,r,Bpack,GAMMALAW);
    flux_set_flux( theFlux , r , n,GAMMALAW,DIVB_CH );
    flux_addto_flux(theFlux,w,grid_NUM_Q(theGrid));
  }

  /*
     if( INCLUDE_VISCOSITY == 1 ){
     double VFlux[NUM_Q];
     double AvgPrim[NUM_Q];
     double Gprim[NUM_Q];
     for( q=0 ; q<NUM_Q ; ++q ){
     AvgPrim[q] = .5*(primL[q]+primR[q]);
     Gprim[q] = .5*(cL->gradp[q]+cR->gradp[q])/r;
     VFlux[q] = 0.0;
     }
     visc_flux( AvgPrim , Gprim , VFlux , r , n );
     for( q=0 ; q<NUM_Q ; ++q ){
     Flux[q] += VFlux[q];
     }
     }
     */
  for( q=0 ; q<NUM_Q ; ++q ){
    cell_add_cons(cL,q,-dt*dA*Flux[q]);
    cell_add_cons(cR,q,dt*dA*Flux[q]);
  }

  cell_add_divB(cL,dA*Bp_face);
  cell_add_divB(cR,-dA*Bp_face);
  cell_add_GradPsi(cL,1,Psi_face/r);
  cell_add_GradPsi(cR,1,-Psi_face/r);

  //free(primL);
  //free(primR);
  //free(Fk);
  free(Uk);
  //free(Flux);
}


