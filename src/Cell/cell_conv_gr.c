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
    g = metric_create(time_global, r, 0, 0);
    a = metric_lapse(g);
    for(i=0; i<3; i++)
        b[i] = metric_shift_u(g, i);
    sqrtg = metric_sqrtgamma(g)/r;

    //Calculate 4-velocity, u0 = u^0, u[i] = u_i
    u0 = 1.0 / sqrt(-metric_g_dd(g,0,0) - 2*metric_dot3_u(g,b,v) - metric_square3_u(g,v));
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

    metric_destroy(g);
}

void cell_cons2prim_gr(double *cons, double *prim, double r, double dV, struct Sim *theSim)
{
    int i,j;
    struct Metric *g;
    double a, b[3], sqrtg;
    double GAMMALAW, rhoh, u0, hmo;
    double rho, v[3], Pp;
    double rhostar, S[3], tau;
    double w, w0, w1; //w = lapse * u^0
    double e, s2, gam; //e = (tau+rhostar)/rhostar, s2 = S_i*S^i/rhostar^2, gam = (GAMMALAW-1)/GAMMALAW
    double c[5];
    double eps = 1.0e-10;
    double CS_FLOOR, CS_CAP, RHO_FLOOR;

    g = metric_create(time_global, r, 0, 0);

    //Metric quantities needed later
    a = metric_lapse(g);
    for(i=0; i<3; i++)
        b[i] = metric_shift_u(g,i);
    sqrtg = metric_sqrtgamma(g)/r;

    //Conserved Quantities
    rhostar = cons[DDD]/dV;
    S[0] = cons[SRR]/dV;
    S[1] = cons[LLL]/dV;
    S[2] = cons[SZZ]/dV;
    tau = cons[TAU]/dV;
    GAMMALAW = sim_GAMMALAW(theSim);

    //Dimensionless conserved variables for w recovery
    e = (tau+rhostar)/rhostar;
    s2 = metric_square3_d(g,S)/(rhostar*rhostar);
    gam = (GAMMALAW-1.0)/GAMMALAW;

    //Coefficients of polynomial to solve for w
    c[0] = -(s2+1)*gam*gam;
    c[1] = 2*gam*e;
    c[2] = gam*gam - e*e + 2*gam*s2;
    c[3] = -2*gam*e;
    c[4] = e*e-s2;
    if(c[4] <= 0)
        printf("e2 <= s2!: e^2 = %lg, s^2 = %lg\n",e*e,s2);

    //Inital guess: previous value of prim
    v[0] = prim[URR];
    v[1] = prim[UPP];
    v[2] = prim[UZZ];
    w0 = a / sqrt(-metric_g_dd(g,0,0) - 2*metric_dot3_u(g,b,v) - metric_square3_u(g,v));
    //Newton-Raphson to find w.
    w1 = w0;
    i = 0;
    do
    {
        w = w1;
        w1 = w - (c[0] + c[1]*w + c[2]*w*w + c[3]*w*w*w + c[4]*w*w*w*w)
                    / (c[1] + 2*c[2]*w + 3*c[3]*w*w + 4*c[4]*w*w*w);
        i++;
    }
    while(fabs(w-w1) > eps && i < 10000);
    if(i == 10000)
    {
        printf("NR failed to converge: w0 = %lg, w = %lg, w1 = %lg\n", w0,w,w1);
        printf("Poly coeffs: c[0]=%lg, c[1]=%lg, c[2]=%lg, c[3]=%lg, c[4]=%lg\n",c[0],c[1],c[2],c[3],c[4]);
    }
    w = w1;

    if(w < 1.0)
    {
        printf("ERROR: w < 1 in cons2prim_gr. (w = %lg)\n", w);
        w = 1.0;
    }
    
    //Prim recovery
    u0 = w/a;
    hmo = w*(e-w)/(w*w-gam);
    if(hmo < -1.0)
        printf("ERROR: h < 0 in cons2prim_gr (h-1 = %lg)\n", hmo);
    else if(hmo < 0.0)
        printf("ERROR?: h < 1 in cons2prim_gr (h-1 = %lg)\n", hmo);

    RHO_FLOOR = sim_RHO_FLOOR(theSim);
    CS_FLOOR = sim_CS_FLOOR(theSim);
    CS_CAP = sim_CS_CAP(theSim);

    rho = rhostar/(sqrtg*w);
    if(rho < RHO_FLOOR)
        rho = RHO_FLOOR;
    Pp = gam * rho * hmo;
    if(Pp < CS_FLOOR*CS_FLOOR*rho*(hmo+1)/GAMMALAW)
        Pp = CS_FLOOR*CS_FLOOR*rho*(hmo+1)/GAMMALAW;
    if(Pp > CS_CAP*CS_CAP*rho*(hmo+1)/GAMMALAW)
        Pp = CS_CAP*CS_CAP*rho*(hmo+1)/GAMMALAW;

    for(i=0; i<3; i++)
    {
        v[i] = 0.0;
        for(j=0; j<3; j++)
            v[i] += metric_gamma_uu(g,i,j) * S[j];
        v[i] /= u0*rhostar*(hmo+1.0);
        v[i] -= b[i];
    }

    prim[RHO] = rho;
    prim[URR] = v[0];
    prim[UPP] = v[1];
    prim[UZZ] = v[2];
    prim[PPP] = Pp;
    
    for(i = sim_NUM_C(theSim); i < sim_NUM_Q(theSim); i++)
        prim[i] = cons[i]/cons[DDD];

    metric_destroy(g);
}

