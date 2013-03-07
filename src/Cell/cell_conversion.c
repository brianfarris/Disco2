#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Grid.h"
#include "../Headers/Face.h"
#include "../Headers/GravMass.h"
#include "../Headers/header.h"

void cell_prim2cons( double * prim , double * cons , double r , double dV ,double GAMMALAW, int runtype){

  double rho = prim[RHO];
  double Pp  = prim[PPP];
  double vr  = prim[URR];
  double vp  = prim[UPP]*r;
  double vz  = prim[UZZ];
  double v2  = vr*vr + vp*vp + vz*vz;
  double rhoe = Pp/(GAMMALAW - 1.);

  cons[DDD] = rho*dV;
  cons[SRR] = rho*vr*dV;
  cons[LLL] = r*rho*vp*dV;
  cons[SZZ] = rho*vz*dV;
  cons[TAU] = ( .5*rho*v2 + rhoe )*dV;

  if (runtype==MHD){
    double Br  = prim[BRR];
    double Bp  = prim[BPP];
    double Bz  = prim[BZZ];
    double B2 = Br*Br + Bp*Bp + Bz*Bz;

    cons[TAU] += .5*B2*dV;
    cons[BRR] = Br*dV/r;
    cons[BPP] = Bp*dV/r;
    cons[BZZ] = Bz*dV;
    cons[PSI] = prim[PSI]*dV;
  }
}

void cell_calc_cons( struct Cell *** theCells,struct Grid *theGrid ){
  double GAMMALAW = grid_GAMMALAW(theGrid);

  int i,j,k;
  for( k=0 ; k<grid_N_z(theGrid) ; ++k ){
    double zp = grid_z_faces(theGrid,k);
    double zm = grid_z_faces(theGrid,k-1);
    double dz = zp-zm;
    for( i=0 ; i<grid_N_r(theGrid) ; ++i ){
      double rp = grid_r_faces(theGrid,i);
      double rm = grid_r_faces(theGrid,i-1);
      double r = .5*(rp+rm);
      for( j=0 ; j<grid_N_p(theGrid,i) ; ++j ){
        struct Cell * c = &(theCells[k][i][j]);
        double dV = .5*(rp*rp - rm*rm)*c->dphi*dz;
        cell_prim2cons( c->prim , c->cons , r , dV,GAMMALAW ,grid_runtype(theGrid));
      }    
    }    
  }
}

void cell_cons2prim( double * cons , double * prim , double r , double dV ,struct Grid * theGrid,int runtype){
  double CS_FLOOR = grid_CS_FLOOR(theGrid);
  double CS_CAP = grid_CS_CAP(theGrid);
  double RHO_FLOOR = grid_RHO_FLOOR(theGrid);
  double VEL_CAP = grid_VEL_CAP(theGrid);
  double GAMMALAW = grid_GAMMALAW(theGrid);
  double rho = cons[DDD]/dV;
  if( rho < RHO_FLOOR ) rho = RHO_FLOOR;
  double Sr  = cons[SRR]/dV;
  double Sp  = cons[LLL]/dV/r;
  double Sz  = cons[SZZ]/dV;
  double E   = cons[TAU]/dV;

  double vr = Sr/rho;
  double vp = Sp/rho;
  double vz = Sz/rho;
  double KE = .5*( Sr*vr + Sp*vp + Sz*vz );
  double v_magnitude = sqrt(2.*KE/rho);

  double B2 = 0.0;
  if (runtype==MHD){
    double Br  = cons[BRR]/dV*r;
    double Bp  = cons[BPP]/dV*r;
    double Bz  = cons[BZZ]/dV;
    B2 = Br*Br+Bp*Bp+Bz*Bz;
    prim[BRR] = Br;
    prim[BPP] = Bp;
    prim[BZZ] = Bz;
    prim[PSI] = cons[PSI]/dV;
  }
  double rhoe = E - KE - .5*B2;
  double Pp = (GAMMALAW-1.)*rhoe;

  if( Pp < CS_FLOOR*CS_FLOOR*rho/GAMMALAW ) {
    Pp = CS_FLOOR*CS_FLOOR*rho/GAMMALAW;
  }
  if ( Pp > CS_CAP*CS_CAP*rho/GAMMALAW) {
    Pp = CS_CAP*CS_CAP*rho/GAMMALAW;
  }
  if ( v_magnitude>VEL_CAP ){
    vr = vr*VEL_CAP/v_magnitude;
    vp = vp*VEL_CAP/v_magnitude;
    vz = vz*VEL_CAP/v_magnitude;
  }
  prim[RHO] = rho;
  prim[PPP] = Pp;
  prim[URR] = vr;
  prim[UPP] = vp/r;
  prim[UZZ] = vz;
}

void cell_calc_prim( struct Cell *** theCells ,struct Grid * theGrid){
  int i,j,k;
  for( k=0 ; k<grid_N_z(theGrid) ; ++k ){
    double zm = grid_z_faces(theGrid,k-1);
    double zp = grid_z_faces(theGrid,k);
    double dz = zp-zm;
    for( i=0 ; i<grid_N_r(theGrid) ; ++i ){
      double rm = grid_r_faces(theGrid,i-1);
      double rp = grid_r_faces(theGrid,i);
      double r = .5*(rp+rm);
      for( j=0 ; j<grid_N_p(theGrid,i) ; ++j ){
        struct Cell * c = cell_single(theCells,i,j,k);
        double dV = .5*(rp*rp-rm*rm)*c->dphi*dz;
        cell_cons2prim( c->cons , c->prim , r , dV ,theGrid,grid_runtype(theGrid));
      }
    }
  }
}



