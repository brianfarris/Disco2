#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Sim.h"
#include "../Headers/Face.h"
#include "../Headers/GravMass.h"
#include "../Headers/header.h"

void cell_prim2cons( double * prim , double * cons , double r , double dV ,struct Sim * theSim){
  double GAMMALAW = sim_GAMMALAW(theSim);
  int runtype = sim_runtype(theSim);

  double rho = prim[RHO];
  double Pp  = prim[PPP];
  double vr  = prim[URR];
  double vp  = prim[UPP]*r;
  double vz  = prim[UZZ];
  double v2  = vr*vr + vp*vp + vz*vz;
  double rhoe = Pp/(GAMMALAW - 1.);

  cons[DDD] = rho*dV;
  cons[SRR] = rho*vr*dV;
  cons[LLL] = r*rho*vp*dV; // note that the conserved variable is angular momentum, not linear momentum
  cons[SZZ] = rho*vz*dV;
  cons[TAU] = ( .5*rho*v2 + rhoe )*dV;

  if (runtype==MHD){
    double Br  = prim[BRR];
    double Bp  = prim[BPP];
    double Bz  = prim[BZZ];
    double B2 = Br*Br + Bp*Bp + Bz*Bz;

    cons[TAU] += .5*B2*dV;
    cons[BRR] = Br*dV;
    cons[BPP] = Bp*dV/r; // note that the conserved variable is Bp/r. This helps us reduce geometric source terms
    cons[BZZ] = Bz*dV;
    cons[PSI] = prim[PSI]*dV;
  }

  int q;
  for( q=sim_NUM_C(theSim) ; q<sim_NUM_Q(theSim) ; ++q ){
    cons[q] = prim[q]*cons[DDD];
  }

}

void cell_calc_cons( struct Cell *** theCells,struct Sim *theSim ){
  double GAMMALAW = sim_GAMMALAW(theSim);

  int i,j,k;
  for( k=0 ; k<sim_N(theSim,Z_DIR) ; ++k ){
    double zp = sim_FacePos(theSim,k,Z_DIR);
    double zm = sim_FacePos(theSim,k-1,Z_DIR);
    double dz = zp-zm;
    for( i=0 ; i<sim_N(theSim,R_DIR) ; ++i ){
      double rp = sim_FacePos(theSim,i,R_DIR);
      double rm = sim_FacePos(theSim,i-1,R_DIR);
      double r = .5*(rp+rm);
      for( j=0 ; j<sim_N_p(theSim,i) ; ++j ){
        struct Cell * c = &(theCells[k][i][j]);
        double dV = .5*(rp*rp - rm*rm)*c->dphi*dz;
        double tau_before = c->cons[TAU];
        cell_prim2cons( c->prim , c->cons , r , dV,theSim);
      }
    }    
  }
}

void cell_cons2prim( double * cons , double * prim , double r , double dV ,struct Sim * theSim){
  int runtype = sim_runtype(theSim);
  double CS_FLOOR = sim_CS_FLOOR(theSim);
  double CS_CAP = sim_CS_CAP(theSim);
  double RHO_FLOOR = sim_RHO_FLOOR(theSim);
  double VEL_CAP = sim_VEL_CAP(theSim);
  double GAMMALAW = sim_GAMMALAW(theSim);
  double rho = cons[DDD]/dV;
  if( rho < RHO_FLOOR ) {
    rho = RHO_FLOOR;
  }
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
    double Br  = cons[BRR]/dV;
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

  int q;
  for( q=sim_NUM_C(theSim) ; q<sim_NUM_Q(theSim) ; ++q ){
    prim[q] = cons[q]/cons[DDD];
  }

}

void cell_calc_prim( struct Cell *** theCells ,struct Sim * theSim){
  int i,j,k;
  for( k=0 ; k<sim_N(theSim,Z_DIR) ; ++k ){
    double zm = sim_FacePos(theSim,k-1,Z_DIR);
    double zp = sim_FacePos(theSim,k,Z_DIR);
    double dz = zp-zm;
    for( i=0 ; i<sim_N(theSim,R_DIR) ; ++i ){
      double rm = sim_FacePos(theSim,i-1,R_DIR);
      double rp = sim_FacePos(theSim,i,R_DIR);
      double r = .5*(rp+rm);
      for( j=0 ; j<sim_N_p(theSim,i) ; ++j ){
        struct Cell * c = cell_single(theCells,i,j,k);
        double dV = .5*(rp*rp-rm*rm)*c->dphi*dz;
        cell_cons2prim( c->cons , c->prim , r , dV ,theSim);
      }
    }
  }
}



