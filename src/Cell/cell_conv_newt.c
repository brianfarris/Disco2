#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Sim.h"
#include "../Headers/Face.h"
#include "../Headers/GravMass.h"
#include "../Headers/header.h"

void cell_prim2cons_newt( double * prim , double * cons , double r , double dV ,struct Sim * theSim){
  double GAMMALAW = sim_GAMMALAW(theSim);

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

  int q;
  for( q=sim_NUM_C(theSim) ; q<sim_NUM_Q(theSim) ; ++q ){
    cons[q] = prim[q]*cons[DDD];
  }

}

void cell_cons2prim_newt( double * cons , double * prim , double r , double dV ,struct Sim * theSim){
  double CS_FLOOR = sim_CS_FLOOR(theSim);
  double CS_CAP = sim_CS_CAP(theSim);
  double RHO_FLOOR = sim_RHO_FLOOR(theSim);
  double VEL_CAP = sim_VEL_CAP(theSim);
  double GAMMALAW = sim_GAMMALAW(theSim);
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

  double rhoe = E - KE;
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

