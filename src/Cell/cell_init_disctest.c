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

// Disc Tests

//Rigidly rotating disc (Flat Space)
void cell_init_disctest_rigid(struct Cell *c, double r, double phi, double z, struct Sim *theSim)
{
    double rho = 1.0;
    double up = sim_InitPar1(theSim);
    double P0 = sim_InitPar2(theSim);
    double GAM = sim_GAMMALAW(theSim);

    double vp = up/sqrt(1.0+r*r*up*up);
    double P = (P0+(GAM-1.0)/GAM*rho)*exp(0.5*GAM*r*r*up*up/(GAM-1.0)) - (GAM-1.0)/GAM*rho;

    if(sim_Background(theSim) == NEWTON)
    {
        vp = up;
        P = P0 + 0.5*rho*r*r*vp*vp;
    }

    c->prim[RHO] = rho;
    c->prim[PPP] = P;
    c->prim[URR] = 0.0;
    c->prim[UPP] = vp;
    c->prim[UZZ] = 0.0;

    if(sim_NUM_C(theSim)<sim_NUM_Q(theSim)) 
    {
        int i;
        for(i=sim_NUM_C(theSim); i<sim_NUM_Q(theSim); i++)
        {
            if(r*cos(phi) < 0)
                c->prim[i] = 0.0;
            else
                c->prim[i] = 1.0;
        }
    }
}

//Keplerian Disc, with optional gaussian bump.
void cell_init_disctest_kepler(struct Cell *c, double r, double phi, double z, struct Sim *theSim)
{
    double rho = 1.0;
    double vr, vp, Pp;
    double B = sim_InitPar1(theSim);
    double mach0 = sim_InitPar2(theSim);
    double r0 = sim_InitPar3(theSim);
    double dr = sim_InitPar4(theSim);
    double GAM = sim_GAMMALAW(theSim);
    double M  = sim_GravM(theSim);
    double nu = fabs(sim_AlphaVisc(theSim));


    rho = 1.0 + B*sqrt(r0/r) + exp(-(r-r0)*(r-r0)/(2*dr*dr));
    Pp = (1.0 + B) * M / (r0 * GAM * mach0*mach0);
    vr = - 3*nu/(2*r) * (1.0+(B*sqrt(r0/r) - 2*r*(r-r0)/(dr*dr)*exp(-(r-r0)*(r-r0)/(2*dr*dr)))/rho);
    vp = sqrt(M/(r*r*r));
    
    if(dr < 0)
    {
        rho = 1.0 + B*sqrt(r0/r);
        vr = -3*nu/(2*r) * (1.0/(1.0+B*sqrt(r0/r)));
    }
   
    c->prim[RHO] = rho;
    c->prim[PPP] = Pp;
    c->prim[URR] = vr;
    c->prim[UPP] = vp;
    c->prim[UZZ] = 0.0;

    if(sim_NUM_C(theSim)<sim_NUM_Q(theSim)) 
    {
        int i;
        for(i=sim_NUM_C(theSim); i<sim_NUM_Q(theSim); i++)
        {
            if(r*cos(phi) < 0)
                c->prim[i] = 0.0;
            else
                c->prim[i] = 1.0;
        }
    }
}

void cell_init_disctest_alphakepler(struct Cell *c, double r, double phi, double z, struct Sim *theSim)
{
    double rho, vr, vp, Pp;
    double rho0 = sim_InitPar1(theSim);
    double P0 = sim_InitPar2(theSim);
    double GAM = sim_GAMMALAW(theSim);
    double M  = sim_GravM(theSim);
    double alpha = fabs(sim_AlphaVisc(theSim));

    rho = rho0;
    Pp = P0;

    double nu = 2*alpha*sqrt(GAM)*Pp/(rho)*sqrt(r*r*r/M);

    vr = -1.5*nu/r;
    vp = exp(-M/(r-2*M)) * sqrt(M/(r*r*r));

    while(1.0-2*M/r-vr*vr/(1-2*M/r)-r*r*vp*vp <= 0.0)
        vr /= 2.0;

    c->prim[RHO] = rho;
    c->prim[PPP] = Pp;
    c->prim[URR] = vr;
    c->prim[UPP] = vp;
    c->prim[UZZ] = 0.0;

    if(sim_NUM_C(theSim)<sim_NUM_Q(theSim)) 
    {
        int i;
        for(i=sim_NUM_C(theSim); i<sim_NUM_Q(theSim); i++)
        {
            if(r*cos(phi) < 0)
                c->prim[i] = 0.0;
            else
                c->prim[i] = 1.0;
        }
    }
}

void cell_init_disctest_spread(struct Cell *c, double r, double phi, double z, struct Sim *theSim)
{
    double rho, vr, vp, P;
    double sig0, T0, vrOut, vrIn, vpOut, vpIn, r0, dr;
    double M  = sim_GravM(theSim);

    sig0 = sim_InitPar1(theSim);
    T0 = sim_InitPar2(theSim);
    r0 = sim_InitPar3(theSim);
    dr = sim_InitPar4(theSim);

    vrOut = 0.0;
    vrIn = -2*M/r / (1+2*M/r);
    vpOut = sqrt(M/(r*r*r));
    vpIn = 0.0;

    //Sigmoid weighting.
    double w = 1.0/(1.0+exp(-(r-r0)/dr));
    
    vr = (1-w)*vrIn + w*vrOut;
    vp = (1-w)*vpIn + w*vpOut;

    rho = sig0;
    P = rho*T0;

    //Cutoff at Vi = 0.5c
    
    vr = vrOut;
    vp = vpOut;
    struct Metric *g = metric_create(time_global, r, phi, z, theSim);
    double a, br, bp, grr, gpp;
    double VMAX = 0.5;
    a = metric_lapse(g);
    br = metric_shift_u(g,0);
    bp = metric_shift_u(g,1);
    grr = metric_gamma_dd(g,0,0);
    gpp = metric_gamma_dd(g,1,1);
    metric_destroy(g);
    
    double VR = sqrt(grr)*(vr+br)/a;
    double VP = sqrt(gpp)*(vp+bp)/a;

    VR = VR >  VMAX ?  VMAX : VR;
    VR = VR < -VMAX ? -VMAX : VR;
    VP = VP >  VMAX ?  VMAX : VP;
    VP = VP < -VMAX ? -VMAX : VP;

    vr = a*VR/sqrt(grr) - br;
    vp = a*VP/sqrt(gpp) - bp;
    
    if(sim_Background(theSim) == GRDISC)
    {
        double W = 1.0/sqrt(1.0-VR*VR-VP*VP);
        double u0 = sqrt(1+2*M/r) * W;
        rho = sig0 / (sqrt(r*r*r*T0/M)/u0);
    }


    c->prim[RHO] = rho;
    c->prim[URR] = vr;
    c->prim[UPP] = vp;
    c->prim[UZZ] = 0.0;
    if(sim_Background(theSim) == GRDISC)
        c->prim[TTT] = T0;
    else
        c->prim[PPP] = P;

    if(sim_NUM_C(theSim)<sim_NUM_Q(theSim)) 
    {
        int i;
        for(i=sim_NUM_C(theSim); i<sim_NUM_Q(theSim); i++)
        {
            if(r*cos(phi) < 0)
                c->prim[i] = 0.0;
            else
                c->prim[i] = 1.0;
        }
    }
}

void cell_init_disctest_calc(struct Cell *c, double r, double phi, double z, struct Sim *theSim)
{
    int test_num = sim_InitPar0(theSim);

    if(test_num == 0)
        cell_init_disctest_rigid(c, r, phi, z, theSim);
    else if(test_num == 1)
        cell_init_disctest_kepler(c, r, phi, z, theSim);
    else if(test_num == 2)
        cell_init_disctest_alphakepler(c, r, phi, z, theSim);
    else if(test_num == 3)
        cell_init_disctest_spread(c, r, phi, z, theSim);
    else
        printf("ERROR: cell_init_disctest given bad option.\n");
}

void cell_single_init_disctest(struct Cell *theCell, struct Sim *theSim,int i,int j,int k)
{
    double rm = sim_FacePos(theSim,i-1,R_DIR);
    double rp = sim_FacePos(theSim,i,R_DIR);
    double r = 0.5*(rm+rp);
    double zm = sim_FacePos(theSim,k-1,Z_DIR);
    double zp = sim_FacePos(theSim,k,Z_DIR);
    double z = 0.5*(zm+zp);
    double t = theCell->tiph-.5*theCell->dphi;

    cell_init_disctest_calc(theCell, r, t, z, theSim);
}

void cell_init_disctest(struct Cell ***theCells,struct Sim *theSim,struct MPIsetup * theMPIsetup)
{
    int i, j, k;
    for (k = 0; k < sim_N(theSim,Z_DIR); k++) 
    {
        double zm = sim_FacePos(theSim,k-1,Z_DIR);
        double zp = sim_FacePos(theSim,k,Z_DIR);
        double z = 0.5*(zm+zp);
        
        for (i = 0; i < sim_N(theSim,R_DIR); i++) 
        {
            double rm = sim_FacePos(theSim,i-1,R_DIR);
            double rp = sim_FacePos(theSim,i,R_DIR);
            double r = 0.5*(rm+rp);

            for (j = 0; j < sim_N_p(theSim,i); j++) 
            {
                double t = theCells[k][i][j].tiph-.5*theCells[k][i][j].dphi;
             
                cell_init_disctest_calc(&(theCells[k][i][j]), r, t, z, theSim);
                
                if(PRINTTOOMUCH)
                {
                    printf("(%d,%d,%d) = (%.12lg, %.12lg, %.12lg): (%.12lg, %.12lg, %.12lg, %.12lg, %.12lg)\n", i,j,k,r,t,z,theCells[k][i][j].prim[RHO],theCells[k][i][j].prim[URR],theCells[k][i][j].prim[UPP],theCells[k][i][j].prim[UZZ],theCells[k][i][j].prim[PPP]);
                }
            }
        }
    }
}
