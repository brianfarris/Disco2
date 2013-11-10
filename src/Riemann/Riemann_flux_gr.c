#define RIEMANN_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Riemann.h"
#include "../Headers/Sim.h"
#include "../Headers/Cell.h"
#include "../Headers/Face.h"
#include "../Headers/Metric.h"
#include "../Headers/header.h"


// ********************************************************************************************
// WE REALLY SHOULD IMPROVE THE COMMENTING OF ALL OF THESE ROUTINES. 
// THERE IS SOME COMPLICATED STUFF HERE. 
// LETS CHOOSE A REFERENCE SUCH AS TORO AND IDENTIFY LINES OF CODE WITH EQUATIONS IN THE BOOK.
// ********************************************************************************************

// this routine is only called by riemann_set_vel.
// It is used to find various L/R quantities. 
void LR_speed_gr(double *prim,double r,int * n,double GAMMALAW,double * p_vn,double * p_cf2,double *Fm,double * p_mn){
  double P   = prim[PPP];
  double rho = prim[RHO];
  double vr  =   prim[URR];
  double vp  = r*prim[UPP];
  double vz  =   prim[UZZ];
  double vn  = vr*n[0] + vp*n[1] + vz*n[2];
  double rhoh = rho + GAMMALAW*P/(GAMMALAW-1);
  double cf2  = GAMMALAW*(P/rhoh);
  double mr = rho*vr;
  double mp = rho*vp;
  double mz = rho*vz;
  double mn = mr*n[0]+mp*n[1]+mz*n[2];

  Fm[0] = rho*vr*vn + P*n[0];
  Fm[1] = rho*vp*vn + P*n[1];
  Fm[2] = rho*vz*vn + P*n[2];

  *p_vn = vn;
  *p_cf2 = cf2;
  *p_mn = mn;
}

// Find velocities needed for the Riemann problem
void riemann_set_vel_gr(struct Riemann * theRiemann,struct Sim * theSim,double r,double GAMMALAW){
  double Sl, Sr, Ss;

  double vnL,cf21,mnL,BnL,B2L;
  double FL[3], FmL[3];
  LR_speed_gr(theRiemann->primL,r,theRiemann->n,GAMMALAW,&vnL,&cf21,FmL,&mnL);

  Sl = vnL - sqrt( cf21 );
  Sr = vnL + sqrt( cf21 );

  double vnR,cf22,mnR,BnR,B2R;
  double FR[3],FmR[3];
  LR_speed_gr(theRiemann->primR,r,theRiemann->n,GAMMALAW,&vnR,&cf22,FmR,&mnR);
 
  if( Sl > vnR - sqrt( cf22 ) ) Sl = vnR - sqrt( cf22 );
  if( Sr < vnR + sqrt( cf22 ) ) Sr = vnR + sqrt( cf22 );
 
  double  mr = ( -Sl*theRiemann->primL[RHO]*theRiemann->primL[URR] + Sr*theRiemann->primR[RHO]*theRiemann->primR[URR] + FmL[0] - FmR[0] )/( Sr - Sl );
  double  mp = ( -Sl*theRiemann->primL[RHO]*theRiemann->primL[UPP]*r + Sr*theRiemann->primR[RHO]*theRiemann->primR[UPP]*r + FmL[1] - FmR[1] )/( Sr - Sl );
  double  mz = ( -Sl*theRiemann->primL[RHO]*theRiemann->primL[UZZ] + Sr*theRiemann->primR[RHO]*theRiemann->primR[UZZ] + FmL[2] - FmR[2] )/( Sr - Sl );
  double rho = ( -Sl*theRiemann->primL[RHO] + Sr*theRiemann->primR[RHO] + mnL - mnR )/( Sr - Sl );

  Ss = (theRiemann->primR[RHO]*vnR*(Sr-vnR)-theRiemann->primL[RHO]*vnL*(Sl-vnL)+theRiemann->primL[PPP]-theRiemann->primR[PPP])
    /(theRiemann->primR[RHO]*(Sr-vnR)-theRiemann->primL[RHO]*(Sl-vnL));

  theRiemann->Sl = Sl;
  theRiemann->Sr = Sr;
  theRiemann->Ss = Ss;

}


//THIS FUNCTION ASSUMES theRiemann->n EQUALS SOME PERMUTATION OF (0,0,1).
//IT WILL NOT WORK FOR ARBITRARY n.
void riemann_set_flux_gr(struct Riemann *theRiemann, struct Sim *theSim, double GAMMALAW, int SetState)
{
    double r = theRiemann->r;
    double *prim;
    double *F;
    
    if (SetState==LEFT)
    {
        prim = theRiemann->primL;
        F = theRiemann->FL;
    }
    else if (SetState==RIGHT)
    {
        prim = theRiemann->primR;
        F = theRiemann->FR;
    } 
    else
    {
        printf("ERROR: riemann_set_flux given unrecognized state.\n");
        exit(0);
    }

    int i,j;
    struct Metric *g;
    double a, b[3], sqrtg;
    double u0, u[3]; //u0 = u^0, u[i] = u_i
    double rho, Pp, v[3];
    double rhoh, vn, bn;

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

    //Calculate 4-velocity
    u0 = sqrt((-1.0 - 2*metric_dot3_u(g, b, v) - metric_square3_u(g,v)) / metric_g_dd(g,0,0));
    for(i=0; i<3; i++)
    {
        u[i] = 0;
        for(j=0; j<3; j++)
            u[i] += metric_gamma_dd(g,i,j) * (v[j]+b[j]);
        u[i] *= u0;
    }

    //Calculate beta & v normal to face.
    vn = v[0]*theRiemann->n[0] + v[1]*theRiemann->n[1] + v[2]*theRiemann->n[2];
    bn = b[0]*theRiemann->n[0] + b[1]*theRiemann->n[1] + b[2]*theRiemann->n[2];
    rhoh = rho+GAMMALAW*Pp/(GAMMALAW-1.);
    
    //Fluxes
    F[DDD] = sqrtg*a*u0 * rho*vn;
    F[SRR] = sqrtg*a*(u0*rhoh * u[0]*vn + Pp*theRiemann->n[0]);
    F[LLL] = sqrtg*a*(u0*rhoh * u[1]*vn + Pp*theRiemann->n[1]);
    F[SZZ] = sqrtg*a*(u0*rhoh * u[2]*vn + Pp*theRiemann->n[2]);
    F[TAU] = sqrtg*(a*u0*(a*u0*rhoh - rho)*vn + Pp*bn);

    //Passive Fluxes
    int q;
    for( q=sim_NUM_C(theSim) ; q<sim_NUM_Q(theSim) ; ++q )
        F[q] = prim[q]*F[DDD];

    metric_destroy(g);
}

