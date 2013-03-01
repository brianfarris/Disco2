#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Grid.h"
#include "../Headers/Riemann.h"
#include "../Headers/GravMass.h"
#include "../Headers/header.h"

void cell_riemann_p( struct Cell * cL , struct Cell * cR, struct Grid * theGrid, double dA , double dt , double r ){
  int NUM_Q = grid_NUM_Q(theGrid);
  double GAMMALAW = grid_GAMMALAW(theGrid);
  double DIVB_CH = grid_DIVB_CH(theGrid);

  struct Riemann * theRiemann = riemann_create(theGrid);

  int q;
  for (q=0;q<NUM_Q;++q){
    riemann_set_primL(theRiemann,q,cell_prims(cL)[q] + cell_gradp(cL)[q]*cell_dphi(cL));
    riemann_set_primR(theRiemann,q,cell_prims(cR)[q] + cell_gradp(cR)[q]*cell_dphi(cR));
  }

  double Sl,Sr,Ss;
  double n[3];
  n[0] = 0.0;   n[1] = 1.0;   n[2] = 0.0;

  double Bpack[6];
  riemann_set_vel(theRiemann,n,r,Bpack,GAMMALAW,DIVB_CH);

  double Bp_face = 0.5*(riemann_prim(theRiemann,LEFT)[BPP]+riemann_prim(theRiemann,RIGHT)[BPP]);
  double Psi_face = 0.5*(riemann_prim(theRiemann,LEFT)[PSI]+riemann_prim(theRiemann,RIGHT)[PSI]);

  if( grid_MOVE_CELLS(theGrid) == C_WRIEMANN ) cell_add_wiph(cL,Ss);
  double w = cell_wiph(cL);
  
  riemann_set_state(theRiemann,w);

  cell_prim2cons( riemann_prim(theRiemann,riemann_state(theRiemann)) , riemann_Uk(theRiemann) , r , 1.0 ,GAMMALAW);
  if((riemann_state(theRiemann)==LEFTSTAR)||(riemann_state(theRiemann)==RIGHTSTAR)){
    riemann_set_Ustar(theRiemann,n,r,Bpack,GAMMALAW,riemann_state(theRiemann));
  }
  riemann_set_flux( theRiemann , r , n,GAMMALAW,DIVB_CH );
  riemann_addto_flux_general(theRiemann,w,grid_NUM_Q(theGrid));


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
    cell_add_cons(cL,q,-dt*dA*riemann_F(theRiemann)[q]);
    cell_add_cons(cR,q,dt*dA*riemann_F(theRiemann)[q]);
  }

  cell_add_divB(cL,dA*Bp_face);
  cell_add_divB(cR,-dA*Bp_face);
  cell_add_GradPsi(cL,1,Psi_face/r);
  cell_add_GradPsi(cR,1,-Psi_face/r);
}


