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

void riemann_visc_flux(struct Riemann *theRiemann, struct Sim *theSim)
{
    int i,j,k,dir;
    double a, sqrtg, du0, hn, cs2, rhoh, visc, height, r;
    double v[3], dv[9], b[3], u[4], du[16], shear[16];
    
    double alpha = sim_AlphaVisc(theSim);
    double GAMMA = sim_GAMMALAW(theSim);
    double M = sim_GravM(theSim);
    int NUMQ = sim_NUM_Q(theSim);
    double *prim = (double *)malloc(NUMQ * sizeof(double));
    double *F = (double *)malloc(NUMQ * sizeof(double));
    struct Metric *g = metric_create(time_global, theRiemann->pos[R_DIR], theRiemann->pos[P_DIR], theRiemann->pos[Z_DIR], theSim);

    r = theRiemann->pos[R_DIR];
    for(i=0; i<NUMQ; i++)
        prim[i] = 0.5*(theRiemann->primL[i] + theRiemann->primR[i]);

    v[0] = prim[URR];
    v[1] = prim[UPP];
    v[2] = prim[UZZ];

    //TODO: FIX THIS.  Doesn't account for misaligned zones in r/z faces.
    dv[0] = 0.5*(cell_gradr(theRiemann->cL, URR) + cell_gradr(theRiemann->cR, URR));
    dv[1] = 0.5*(cell_gradr(theRiemann->cL, UPP) + cell_gradr(theRiemann->cR, UPP));
    dv[2] = 0.5*(cell_gradr(theRiemann->cL, UZZ) + cell_gradr(theRiemann->cR, UZZ));
    dv[3] = 0.5*(cell_gradp(theRiemann->cL, URR) + cell_gradp(theRiemann->cR, URR));
    dv[4] = 0.5*(cell_gradp(theRiemann->cL, UPP) + cell_gradp(theRiemann->cR, UPP));
    dv[5] = 0.5*(cell_gradp(theRiemann->cL, UZZ) + cell_gradp(theRiemann->cR, UZZ));
    dv[6] = 0.5*(cell_gradz(theRiemann->cL, URR) + cell_gradz(theRiemann->cR, URR));
    dv[7] = 0.5*(cell_gradz(theRiemann->cL, UPP) + cell_gradz(theRiemann->cR, UPP));
    dv[8] = 0.5*(cell_gradz(theRiemann->cL, UZZ) + cell_gradz(theRiemann->cR, UZZ));

    a = metric_lapse(g);
    for(i=0; i<3; i++)
        b[i] = metric_shift_u(g, i);
    sqrtg = metric_sqrtgamma(g)/r;

    u[0] = 1.0/sqrt(-metric_g_dd(g,0,0) - 2*metric_dot3_u(g,b,v) - metric_square3_u(g,v));
    for(i=0; i<3; i++)
        u[i+1] = u[0]*v[i];

    du[0] = 0; du[1] = 0; du[2] = 0; du[3] = 0;
    for(i=1; i<4; i++)
    {
        du0 = metric_dg_dd(g,i,0,0);
        for(j=1; j<3; j++)
        {
            du0 += 2*v[j-1]*metric_dg_dd(g,i,0,j) + 2*metric_g_dd(g,0,j)*dv[3*(i-1)+j-1];
            for(k=1; k<4; k++)
                du0 += v[j-1]*v[k-1]*metric_dg_dd(g,i,j,k) + 2*metric_g_dd(g,j,k)*v[j-1]*dv[3*(i-1)+j-1];
        }
        du[4*i] = 0.5*u[0]*u[0]*u[0]*du0;
    }

    for(i=1; i<4; i++)
        for(j=1; j<4; j++)
            du[4*i+j] = v[j-1]*du[4*i] + u[0]*dv[3*(i-1)+j-1];

    metric_shear_uu(g, u, du, shear);

    //TODO: CHECK THIS!  Probably not relativistic and/or isothermal.
    rhoh = prim[RHO] + GAMMA*prim[PPP]/(GAMMA-1.0);
    cs2 = GAMMA*prim[PPP] / rhoh;
    height = sqrt(prim[PPP]*r*r*r*(1-3*M/r)/(rhoh*M));
    if(alpha > 0)
        visc = alpha * sqrt(cs2) * height;
    else
        visc = -alpha;

    for(i=0; i<3; i++)
        if(theRiemann->n[i] == 1)
            dir = i;

    if(dir == 1)
        hn = theRiemann->pos[R_DIR];
    else
        hn = 1.0;

    for(i=0; i<NUMQ; i++)
        F[i] = 0.0;

    for(i=0; i<4; i++)
    {
        F[SRR] += metric_g_dd(g,i,1)*shear[4*(dir+1)+i];
        F[LLL] += metric_g_dd(g,i,2)*shear[4*(dir+1)+i];
        F[SZZ] += metric_g_dd(g,i,3)*shear[4*(dir+1)+i];
    }
    F[SRR] *= -a*sqrtg*hn * visc;
    F[LLL] *= -a*sqrtg*hn * visc;
    F[SZZ] *= -a*sqrtg*hn * visc;
    F[TAU] = -a*a*sqrtg*hn * visc * shear[dir+1];

    theRiemann->F[SRR] += F[SRR];
    theRiemann->F[LLL] += F[LLL];
    theRiemann->F[SZZ] += F[SZZ];
    theRiemann->F[TAU] += F[TAU];

    free(prim);
    free(F);
    metric_destroy(g);
}

