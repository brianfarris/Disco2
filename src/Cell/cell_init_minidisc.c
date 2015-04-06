#define CELL_PRIVATE_DEFS
#define EOS_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/EOS.h"
#include "../Headers/Face.h"
#include "../Headers/GravMass.h"
#include "../Headers/header.h"
#include "../Headers/Metric.h"
#include "../Headers/Sim.h"

//Minidisc initial condition files.

void cell_init_minidisc_spread(struct Cell *c, double r, double phi, double z,
                                struct Sim *theSim)
{
    int i;
    double M = sim_GravM(theSim);
    double alpha = sim_AlphaVisc(theSim);

    double rho, T, vr, vp, q;
    double R0, dr, rho1, rho0, T0;

    R0 = sim_InitPar1(theSim);   //Location of disc
    dr = sim_InitPar2(theSim);  //Extent of disc from R0
    rho1 = sim_InitPar3(theSim); //Peak density at R0
    rho0 = sim_InitPar4(theSim); //Density of atmosphere
    T0 = sim_InitPar5(theSim);   //Temperature of atmosphere

    struct Metric *g = metric_create(time_global, r, phi, z, theSim);
    double a, b[3];
    a = metric_lapse(g);
    b[0] = metric_shift_u(g,0);
    b[1] = metric_shift_u(g,1);
    b[2] = metric_shift_u(g,2);
    metric_destroy(g);

    rho = rho0;
    T = T0;
    vr = -b[0];
    vp = -b[1];
    q = 0.0;

    if( fabs(r-R0) < 5*dr)
    {
        double fac = exp(-(r-R0)*(r-R0)/(2*dr*dr));
        //First get everything in CGS.
        rho += fac * rho1;
        vr = 0.0;
        vp += fac*sqrt(M/(r*r*r));
        q = 1;
    }
/*
    if(sim_Background(theSim) != GRDISC)
    {
        double tp[5];
        tp[RHO] = rho0;
        tp[TTT] = T0;
        tp[URR] = vr;
        tp[UPP] = vp;
        tp[UZZ] = -b[2];
        double eps = eos_eps(tp, theSim);
        double P = eos_ppp(tp, theSim);
        double H = sqrt(r*r*r*P/(M*(rho+rho*eps+P))) * a;
        rho = rho0 * H;
        T = P * H;
    }
*/
    c->prim[RHO] = rho;
    c->prim[TTT] = T;
    c->prim[URR] = vr;
    c->prim[UPP] = vp;
    c->prim[UZZ] = -b[2];

    for(i=sim_NUM_C(theSim); i<sim_NUM_Q(theSim); i++)
        c->prim[i] = q;
}
void cell_init_minidisc_shakura(struct Cell *c, double r, double phi, double z,
                            struct Sim *theSim)
{
    int i;
    double M = sim_GravM(theSim);
    double alpha = sim_AlphaVisc(theSim);

    double rho, T, vr, vp, q;
    double R0, R1, mdot, rho0, T0;

    R0 = sim_InitPar1(theSim);   //Inner Edge of Disc
    R1 = sim_InitPar2(theSim);  //Outer Edge of Disc
    mdot = sim_InitPar3(theSim); //Accretion rate in solar masses per second
    rho0 = sim_InitPar4(theSim); //Density of atmosphere
    T0 = sim_InitPar5(theSim);   //Temperature of atmosphere

    struct Metric *g = metric_create(time_global, r, phi, z, theSim);
    double a, b[3];
    a = metric_lapse(g);
    b[0] = metric_shift_u(g,0);
    b[1] = metric_shift_u(g,1);
    b[2] = metric_shift_u(g,2);
    metric_destroy(g);

    rho = rho0;
    T = T0;
    vr = -b[0];
    vp = -b[1];
    q = 0.0;

    //if(r > 10*M)
    //{
        //First get everything in CGS.
        double omk = sqrt(M/(r*r*r)) * eos_c/eos_r_scale;
        double kappa = 0.2;  //Thompson opacity in cgs
        mdot *= 1.9884e33;  

        double PH = mdot/(3*M_PI*alpha) * omk;
        double qdot = 3*omk*omk*mdot/(3*M_PI);
        double rhoH = pow(0.75 * eos_sb*eos_mp*eos_mp*eos_mp*eos_mp/kappa
                             * PH*PH*PH*PH / qdot, 0.2);
        double vrcgs = -mdot / (2*M_PI*r*rhoH);

        //cavity scale
        double scale = exp(-pow(r/R0,-1) - pow(r/R1,10));
        //double scale = exp(-10*(R0/r + R1/r));
        
        //Now assign values in code units
        rho = scale * rhoH / (eos_rho_scale*eos_r_scale) + (1-scale)*rho0;
        T = scale * PH / (eos_rho_scale*eos_c*eos_c*eos_r_scale) + (1-scale)*T0;
        vr = vrcgs / eos_c;
        vp = omk / (eos_c/eos_r_scale);
        q += scale * 1.0;
    //}

    double rel_weight = exp(-pow(r/(6*M),-5));
    vr = rel_weight * vr + (1-rel_weight) * (-b[0]);
    if(vr < (-2*M/r) / (1+2*M/r))
        vr = (-2*M/r) / (1+2*M/r);
    vp = rel_weight * vp + (1-rel_weight) * 1.0*sqrt(M/(r*r*(r+2*M)));
/*
    if(sim_Background(theSim) != GRDISC)
    {
        double tp[5];
        tp[RHO] = rho0;
        tp[TTT] = T0;
        tp[URR] = vr;
        tp[UPP] = vp;
        tp[UZZ] = -b[2];
        double eps = eos_eps(tp, theSim);
        double P = eos_ppp(tp, theSim);
        double H = sqrt(r*r*r*P/(M*(rho+rho*eps+P))) * a;
        rho = rho0 * H;
        T = P * H;
    }
*/
    c->prim[RHO] = rho;
    c->prim[TTT] = T;
    c->prim[URR] = vr;
    c->prim[UPP] = vp;
    c->prim[UZZ] = -b[2];

    for(i=sim_NUM_C(theSim); i<sim_NUM_Q(theSim); i++)
        c->prim[i] = q;
}

void cell_init_minidisc_calc(struct Cell *c, double r, double phi, double z,
                            struct Sim *theSim)
{
    int opt = sim_InitPar0(theSim);

    if(opt == 0)
        cell_init_minidisc_spread(c, r, phi, z, theSim);
    else if(opt == 1)
        cell_init_minidisc_shakura(c, r, phi, z, theSim);
    else
        printf("ERROR: cell_init_minidisc given bad option.\n");
}

void cell_single_init_minidisc(struct Cell *theCell, struct Sim *theSim,
                            int i, int j, int k)
{
    double rm = sim_FacePos(theSim,i-1,R_DIR);
    double rp = sim_FacePos(theSim,i,R_DIR);
    double r = 0.5*(rm+rp);
    double zm = sim_FacePos(theSim,k-1,Z_DIR);
    double zp = sim_FacePos(theSim,k,Z_DIR);
    double z = 0.5*(zm+zp);
    double t = theCell->tiph-.5*theCell->dphi;

    cell_init_minidisc_calc(theCell, r, t, z, theSim);
}

void cell_init_minidisc(struct Cell ***theCells,struct Sim *theSim,
                    struct MPIsetup *theMPIsetup)
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

                cell_init_minidisc_calc(&(theCells[k][i][j]), r, t, z, theSim);

                if(PRINTTOOMUCH)
                {
                    printf("(%d,%d,%d) = (%.12lg, %.12lg, %.12lg): (%.12lg, %.12lg, %.12lg, %.12lg, %.12lg)\n",
                          i,j,k,r,t,z,theCells[k][i][j].prim[RHO], theCells[k][i][j].prim[URR],
                          theCells[k][i][j].prim[UPP], theCells[k][i][j].prim[UZZ],
                          theCells[k][i][j].prim[PPP]);
                }
            }
        }
    }
}
