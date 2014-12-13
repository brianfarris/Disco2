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

void cell_prim2cons_grdisc(double *prim, double *cons, double *pos, 
                            double dV, struct Sim *theSim)
{
    struct Metric *g;
    double a, b[3], sqrtg;
    double u0, u_d[4], U[4];
    double rho, T, v[3];
    double Pp, eps, rhoh, H, M;
    double r, phi, z;
    int i,j;

    //Get hydro primitives
    rho = prim[RHO];
    v[0] = prim[URR];
    v[1] = prim[UPP];
    v[2] = prim[UZZ];

    //Get ADM metric values
    g = metric_create(time_global, r, phi, z, theSim);
    a = metric_lapse(g);
    for(i=0; i<3; i++)
        b[i] = metric_shift_u(g, i);
    sqrtg = metric_sqrtgamma(g) / r;
    for(i=0; i<4; i++)
        U[i] = metric_frame_U_u(g, i, theSim);
    M = sim_GravM(theSim);

    //Calculate 4-velocity
    u0 = 1.0 / sqrt(-metric_g_dd(g,0,0) - 2*metric_dot3_u(g,b,v)
                    - metric_square3_u(g,v));
    u_d[0] = metric_g_dd(g,0,0)*u0 + metric_dot3_u(g,b,v)*u0;
    for(i=1; i<4; i++)
    {
        u_d[i] = 0.0;
        for(j=0; j<3; j++)
            u_d[i] += metric_gamma_dd(g,i-1,j) * (v[j] + b[j]);
        u_d[i] *= u0;
    }

    Pp = eos_ppp(prim, theSim);
    eps = eos_eps(prim, theSim);
    rhoh = rho*(1.0+eps) + Pp;

    H = sqrt(r*r*r*Pp / (rhoh*M)) / u0;

    cons[DDD] = a*sqrtg*u0 * rho * H * dV;
    cons[SRR] = a*sqrtg*rhoh*u0 * u_d[1] * H * dV;
    cons[LLL] = a*sqrtg*rhoh*u0 * u_d[2] * H * dV;
    cons[SZZ] = a*sqrtg*rhoh*u0 * u_d[3] * H * dV;
    cons[TAU] = a*sqrtg*(-rhoh*u0*(U[0]*u_d[0]+U[1]*u_d[1]+U[2]*u_d[2]
                        +U[3]*u_d[3]) - rho*u0 - U[0]*Pp) * H * dV;

    for(i=sim_NUM_C(theSim); i<sim_NUM_Q(theSim); i++)
        cons[i] = prim[i] * cons[DDD];

    metric_destroy(g);
}

//TODO: WRITE THIS
void cell_cons2prim_grdisc(double *cons, double *prim, double *pos, double dV, struct Sim *theSim)
{
}

