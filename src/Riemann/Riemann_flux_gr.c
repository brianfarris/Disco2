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
// TODO: remove GAMMALAW
void LR_speed_gr(double *prim, double r, int *n ,double GAMMALAW, double *p_vn, double *p_cf2, double *Fm, double *p_mn)
{
    double P   = prim[PPP];
    double rho = prim[RHO];
    double vr  = prim[URR];
    double vp  = prim[UPP];
    double vz  = prim[UZZ];
    double vn  = vr*n[0] + vp*n[1] + vz*n[2];
    double cf2  = GAMMALAW*P/(rho + GAMMALAW*P/(GAMMALAW-1.0));

    *p_vn = vn;
    *p_cf2 = cf2;
    
    //TODO: These are only used for HLLC and are WRONG.
    *p_mn = 0;
    Fm[0] = 0;
    Fm[1] = 0;
    Fm[2] = 0;
}

// Find velocities needed for the Riemann problem
void riemann_set_vel_gr(struct Riemann *theRiemann, struct Sim *theSim, double r, double GAMMALAW)
{
    int i, dir;
    double Sl, Sr, Sl1, Sr1, Sl2, Sr2;
    double a, b[3], bn, sig, gam, w, v[3], dv;
    struct Metric *g;

    g = metric_create(time_global, theRiemann->pos[R_DIR], theRiemann->pos[P_DIR], theRiemann->pos[Z_DIR], theSim);

    dir = -1;
    for(i=0; i<3; i++)
    {
        b[i] = metric_shift_u(g,i);
        if(theRiemann->n[i] != 0)
            dir = i;
    }
    a = metric_lapse(g);
    bn = b[dir];
    gam = metric_gamma_uu(g, dir, dir);

    double vnL, cf21, mnL;
    double FmL[3];
    LR_speed_gr(theRiemann->primL, r, theRiemann->n, GAMMALAW, &vnL, &cf21, FmL, &mnL);

    v[0] = theRiemann->primL[URR];
    v[1] = theRiemann->primL[UPP];
    v[2] = theRiemann->primL[UZZ];
    w = a / sqrt(-metric_g_dd(g,0,0)-2.0*metric_dot3_u(g,b,v)-metric_square3_u(g,v));
    sig = cf21/(w*w*(1.0-cf21));
    dv = sqrt(sig*(1.0+sig)*a*a*gam - sig*(vnL+bn)*(vnL+bn));
    
    Sl1 = (vnL - sig*bn - dv) / (1.0+sig);
    Sr1 = (vnL - sig*bn + dv) / (1.0+sig);

    //TODO: Remove: Lax-Friedrichs
//    Sl1 = -a*sqrt(gam) - bn;
//    Sr1 = a*sqrt(gam) - bn;
//    Sl1 = -1.0;
//    Sr1 = 1.0;

    double vnR,cf22,mnR;
    double FmR[3];
    LR_speed_gr(theRiemann->primR, r, theRiemann->n, GAMMALAW, &vnR, &cf22, FmR, &mnR);

    v[0] = theRiemann->primR[URR];
    v[1] = theRiemann->primR[UPP];
    v[2] = theRiemann->primR[UZZ];
    w = a / sqrt(-metric_g_dd(g,0,0)-2.0*metric_dot3_u(g,b,v)-metric_square3_u(g,v));
    sig = cf22/(w*w*(1.0-cf22));
    dv = sqrt(sig*(1.0+sig)*a*a*gam - sig*(vnL+bn)*(vnL+bn));

    Sl2 = (vnR - sig*bn - dv) / (1.0+sig);
    Sr2 = (vnR - sig*bn + dv) / (1.0+sig);

    //Lax-Friedrichs
//    Sl2 = -a*sqrt(gam) - bn;
//    Sr2 = a*sqrt(gam) - bn;
//    Sl2 = -1.0;
//    Sr2 = 1.0;

    if(Sl1 > Sl2)
        Sl = Sl2;
    else
        Sl = Sl1;
    if(Sr1 < Sr2)
        Sr = Sr2;
    else
        Sr = Sr1;

    //Fluxes in orthonormal basis
    if(dir == PDIRECTION)
    {
        Sl *= r;
        Sr *= r;
    }

  theRiemann->Sl = Sl;
  theRiemann->Sr = Sr;

  //TODO: This is only used for HLLC and is WRONG.
  theRiemann->Ss = (Sl+Sr)/2;

  metric_destroy(g);
}

//THIS FUNCTION ASSUMES theRiemann->n EQUALS SOME PERMUTATION OF (0,0,1).
//IT WILL NOT WORK FOR ARBITRARY n.
void riemann_set_flux_gr(struct Riemann *theRiemann, struct Sim *theSim, double GAMMALAW, int SetState)
{
    double r = theRiemann->pos[R_DIR];
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
    double a, b[3], sqrtg, U[4];
    double u0, u[4]; //u0 = u^0, u[i] = u_i
    double rho, Pp, v[3];
    double rhoh, vn, bn, hn, Un;

    //Get hydro primitives
    rho = prim[RHO];
    Pp  = prim[PPP];
    v[0]  = prim[URR];
    v[1]  = prim[UPP];
    v[2]  = prim[UZZ];

    //Get needed metric values
    g = metric_create(time_global, theRiemann->pos[R_DIR], theRiemann->pos[P_DIR], theRiemann->pos[Z_DIR], theSim);
    a = metric_lapse(g);
    for(i=0; i<3; i++)
        b[i] = metric_shift_u(g, i);
    sqrtg = metric_sqrtgamma(g)/r;
    for(i=0; i<4; i++)
        U[i] = metric_frame_U_u(g,i,theSim);

    //Check if interpolated velocity is superluminal
    if(-metric_g_dd(g,0,0) - 2*metric_dot3_u(g,b,v) - metric_square3_u(g,v) < 0)
    {
        printf("ERROR: Velocity too high in flux. r=%.12g, vr=%.12g, vp=%.12g, vz=%.12g\n", r, v[0], v[1], v[2]);

        //If velocity is superluminal, reduce to Lorentz factor 5, keeping
        //direction same in rest frame.
        
        double MAXW = 5;
        double V[3], V2, corr;
        //Calculate Eulerian velocity
        for(i=0; i<3; i++)
            V[i] = (v[i]+b[i])/a;
        V2 = metric_square3_u(g, V);
        //Correction factor.
        corr = sqrt((MAXW*MAXW+1.0)/(MAXW*MAXW*V2));
        //Reset velocity
        for(i=0; i<3; i++)
            v[i] = corr*v[i] - (1.0-corr)*b[i];

        printf("   fix: badV2 = %.12g corr = %.12g, newV2 = %.12g\n", 
                V2, corr, 1.0-1.0/(MAXW*MAXW));
    }
    //TODO: Remove this when the fix above is tested.
    if(-metric_g_dd(g,0,0) - 2*metric_dot3_u(g,b,v) - metric_square3_u(g,v) < 0)
    {
        printf("AGAIN.  r = %.12g\n", r);
    }

    //Calculate 4-velocity
    u0 = 1.0 / sqrt(-metric_g_dd(g,0,0) - 2*metric_dot3_u(g,b,v) - metric_square3_u(g,v));
    u[0] = metric_g_dd(g,0,0) * u0 + u0*metric_dot3_u(g,b,v);
    for(i=1; i<4; i++)
    {
        u[i] = 0;
        for(j=0; j<3; j++)
            u[i] += metric_gamma_dd(g,i-1,j) * (v[j]+b[j]);
        u[i] *= u0;
    }

    //Calculate beta & v normal to face.
    vn = v[0]*theRiemann->n[0] + v[1]*theRiemann->n[1] + v[2]*theRiemann->n[2];
    bn = b[0]*theRiemann->n[0] + b[1]*theRiemann->n[1] + b[2]*theRiemann->n[2];
    Un = U[1]*theRiemann->n[0] + U[2]*theRiemann->n[1] + U[3]*theRiemann->n[2];
    if(theRiemann->n[1] == 1)
        hn = r;
    else
        hn = 1.0;
   
    rhoh = rho + GAMMALAW/(GAMMALAW-1.0)*Pp;
    
    //Fluxes
    F[DDD] = hn*sqrtg*a*u0 * rho*vn;
    F[SRR] = hn*sqrtg*a*(u0*rhoh * u[1]*vn + Pp*theRiemann->n[0]);
    F[LLL] = hn*sqrtg*a*(u0*rhoh * u[2]*vn + Pp*theRiemann->n[1]);
    F[SZZ] = hn*sqrtg*a*(u0*rhoh * u[3]*vn + Pp*theRiemann->n[2]);
    F[TAU] = hn*sqrtg*a*(u0*(-rhoh*(U[0]*u[0]+U[1]*u[1]+U[2]*u[2]+U[3]*u[3]) - rho)*vn - Un*Pp);
    //F[TAU] = hn*sqrtg*(a*u0*(a*u0*rhoh - rho)*vn + Pp*bn);

    //HLL Viscous Flux
    /*
    double *Fvisc = (double *) malloc(sim_NUM_Q(theSim) * sizeof(double));
    riemann_visc_flux_LR(theRiemann, theSim, SetState, Fvisc);
    F[DDD] += Fvisc[DDD];
    F[SRR] += Fvisc[SRR];
    F[LLL] += Fvisc[LLL];
    F[SZZ] += Fvisc[SZZ];
    F[TAU] += Fvisc[TAU];
    free(Fvisc);
    */

    //Passive Fluxes
    int q;
    for( q=sim_NUM_C(theSim) ; q<sim_NUM_Q(theSim) ; ++q )
        F[q] = prim[q]*F[DDD];

    metric_destroy(g);
}

