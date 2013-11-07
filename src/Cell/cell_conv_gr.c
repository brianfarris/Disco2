#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Sim.h"
#include "../Headers/Face.h"
#include "../Headers/GravMass.h"
#include "../Headers/Metric.h"
#include "../Headers/header.h"

void cell_prim2cons_gr( double * prim , double * cons , double r , double dV ,struct Sim * theSim)
{
    int i,j;

    struct Metric *g;
    double a, b[3], sqrtg;
    double u0, u[3];
    double rho, Pp, v[3];
    double GAMMALAW, rhoh;

    //Get hydro primitives
    rho = prim[RHO];
    Pp  = prim[PPP];
    v[0]  = prim[URR];
    v[1]  = prim[UPP];
    v[2]  = prim[UZZ];

    //Get needed metric values
    g = metric_create(theSim, time_global, r, 0, 0);
    a = metric_lapse(g);
    for(i=0; i<3; i++)
        b[i] = metric_shift_u(g, i);
    sqrtg = metric_sqrtgamma(g);

    //Calculate 4-velocity, u0 = u^0, u[i] = u_i
    u0 = (-1.0 - 2*metric_dot3_u(g, b, v) - metric_square3_u(g,v)) / metric_g_dd(g,0,0);
    for(i=0; i<3; i++)
    {
        u[i] = 0;
        for(j=0; j<3; j++)
            u[i] += metric_gamma_dd(g,i,j) * (v[j]+b[j]);
        u[i] *= u0;
    }
    
    GAMMALAW = sim_GAMMALAW(theSim);
    rhoh = rho + GAMMALAW*Pp/(GAMMALAW - 1.);

    cons[DDD] = a*sqrtg*u0 * rho * dV;
    cons[SRR] = a*sqrtg*rhoh*u0 * u[0] * dV;
    cons[LLL] = a*sqrtg*rhoh*u0 * u[1] * dV;
    cons[SZZ] = a*sqrtg*rhoh*u0 * u[2] * dV;
    cons[TAU] = sqrtg * (a*u0*(a*u0*rhoh - rho) - Pp) * dV;

    int q;
    for( q=sim_NUM_C(theSim) ; q<sim_NUM_Q(theSim) ; ++q ){
        cons[q] = prim[q]*cons[DDD];
    }
}

void cell_cons2prim_gr( double * cons , double * prim , double r , double dV ,struct Sim * theSim){
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

