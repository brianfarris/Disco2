#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/EOS.h"
#include "../Headers/Face.h"
#include "../Headers/GravMass.h"
#include "../Headers/header.h"
#include "../Headers/Sim.h"

//ADAF2: Advection-Dominated Accretion Flow for GRDISC

void cell_init_adaf2_nocool(struct Cell *c, double r, double phi, double z,
                            struct Sim *theSim)
{
    int i;
    double M = sim_GravM(theSim);
    double alpha = sim_AlphaVisc(theSim);

    double rho, vr, vp, T;
    double rho0, T0, vr0, r0, u0;
    double rho_power, T_power;
    double vMax = 0.5;

    r0 = sim_InitPar1(theSim);
    rho0 = sim_InitPar2(theSim);
    T0 = sim_InitPar3(theSim);
    vr0 = sim_InitPar4(theSim);
    rho_power = sim_InitPar5(theSim);
    T_power = sim_InitPar6(theSim);

    rho = rho0 * pow(r/r0, rho_power);
    T = T0 * pow(r/r0, T_power);
    vr = vr0 * alpha * pow(r/M, -0.5);
    vp = pow(r/M, -1.5) / M;

    if(sim_Metric(theSim) == SCHWARZSCHILD_SC)
    {
        if(vr > vMax * (1-2*M/r))
            vr = vMax * (1-2*M/r);
        if(vr < -vMax * (1-2*M/r))
            vr = -vMax * (1-2*M/r);
        if(vp > vMax * sqrt(1-2*M/r)/r)
            vp = vMax * sqrt(1-2*M/r)/r;
        if(vp < -vMax * sqrt(1-2*M/r)/r)
            vp = -vMax * sqrt(1-2*M/r)/r;
        u0 = 1.0/sqrt(1-2*M/r-vr*vr/(1-2*M/r)-r*r*vp*vp);
    }
    else if(sim_Metric(theSim) == SCHWARZSCHILD_KS)
    {
        /*
        double fac = 1.0/(1.0 + 2*M*vr/(r-2*M));
        vr *= fac;
        vp *= fac;
        */

        if(vr > (vMax-2*M/r)/(1+2*M/r))
            vr = (vMax-2*M/r)/(1+2*M/r);
        if(vr < (-vMax-2*M/r)/(1+2*M/r))
            vr = (-vMax-2*M/r)/(1+2*M/r);
        if(vp > vMax/(sqrt(1+2*M/r)*r))
            vp = vMax/(sqrt(1+2*M/r)*r);
        if(vp < -vMax/(sqrt(1+2*M/r)*r))
            vp = -vMax/(sqrt(1+2*M/r)*r);
        u0 = 1.0/sqrt(1-2*M/r-4*M/r*vr-(1+2*M/r)*vr*vr-r*r*vp*vp);
    }
    else if(sim_Metric(theSim) == KERR_KS)
    {
        if(vr > (vMax-2*M/r)/(1+2*M/r))
            vr = (vMax-2*M/r)/(1+2*M/r);
        if(vr < (-vMax-2*M/r)/(1+2*M/r))
            vr = (-vMax-2*M/r)/(1+2*M/r);
        if(vp > vMax/(sqrt(1+2*M/r)*r))
            vp = vMax/(sqrt(1+2*M/r)*r);
        if(vp < -vMax/(sqrt(1+2*M/r)*r))
            vp = -vMax/(sqrt(1+2*M/r)*r);
        u0 = 1.0/sqrt(1-2*M/r-4*M/r*vr-(1+2*M/r)*vr*vr-r*r*vp*vp);
    }
    else
        u0 = 1.0;

    if(sim_Background(theSim) != GRDISC)
    {
        double tp[5];
        tp[RHO] = rho;
        tp[TTT] = T;
        tp[URR] = vr;
        tp[UPP] = vp;
        tp[UZZ] = 0.0;
        double eps = eos_eps(tp, theSim);
        double P = eos_ppp(tp, theSim);
        double H = sqrt(r*r*r*P/(M*(rho+rho*eps+P))) / u0;
        rho = rho * H;
        T = P * H;
    }

    c->prim[RHO] = rho;
    c->prim[TTT] = T;
    c->prim[URR] = vr;
    c->prim[UPP] = vp;
    c->prim[UZZ] = 0.0;

    for(i=sim_NUM_C(theSim); i<sim_NUM_Q(theSim); i++)
    {
        if(r*cos(phi) < 0)
            c->prim[i] = 0.0;
        else
            c->prim[i] = 1.0;
    }
}

void cell_init_adaf2_calc(struct Cell *c, double r, double phi, double z,
                            struct Sim *theSim)
{
    int disc_num = sim_InitPar0(theSim);

    if(disc_num == 0)
        cell_init_adaf2_nocool(c, r, phi, z, theSim);
    else
        printf("ERROR: cell_init_accdisc given bad option.\n");
}

void cell_single_init_adaf2(struct Cell *theCell, struct Sim *theSim,
                            int i, int j, int k)
{
    double rm = sim_FacePos(theSim,i-1,R_DIR);
    double rp = sim_FacePos(theSim,i,R_DIR);
    double r = 0.5*(rm+rp);
    double zm = sim_FacePos(theSim,k-1,Z_DIR);
    double zp = sim_FacePos(theSim,k,Z_DIR);
    double z = 0.5*(zm+zp);
    double t = theCell->tiph-.5*theCell->dphi;

    cell_init_adaf2_calc(theCell, r, t, z, theSim);
}

void cell_init_adaf2(struct Cell ***theCells,struct Sim *theSim,
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

                cell_init_adaf2_calc(&(theCells[k][i][j]), r, t, z, theSim);

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
