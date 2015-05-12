#define CELL_PRIVATE_DEFS
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

//Kicked Discs

void cell_init_kick_isothermal(struct Cell *c, double r, double phi, double z,
                                struct Sim *theSim)
{
    int i;
    double M = sim_GravM(theSim);
    double al = sim_AlphaVisc(theSim);
    double GAM = sim_GAMMALAW(theSim);
    double VMAX = 0.5;

    double rho0, T, r0, VR, VP, rho, P, vr, vp, q, kick;

    r0 = sim_InitPar1(theSim);
    rho0 = sim_InitPar2(theSim);
    T = sim_InitPar3(theSim);
    kick = sim_InitPar4(theSim);

    struct Metric *g = metric_create(time_global, r, phi, z, theSim);
    double a, br, bp, grr, gpp;
    a = metric_lapse(g);
    br = metric_shift_u(g,0);
    bp = metric_shift_u(g,1);
    grr = metric_gamma_dd(g,0,0);
    gpp = metric_gamma_dd(g,1,1);
    metric_destroy(g);

    rho = rho0 * pow(r/r0, -1.5);
    P = T * rho;
    vr = -1.5 * al * sqrt(GAM) * T * sqrt(r/M);
    vp = sqrt(M/(r*r*r));
    q = r>r0 ? 0.0 : 1.0;

    VR = sqrt(grr)*(vr+br)/a - kick*sin(phi);
    VP = sqrt(gpp)*(vp+bp)/a - kick*cos(phi)/r;

    VR = VR >  VMAX ?  VMAX : VR;
    VR = VR < -VMAX ? -VMAX : VR;
    VP = VP >  VMAX ?  VMAX : VP;
    VP = VP < -VMAX ? -VMAX : VP;

    vr = a*VR/sqrt(grr) - br;
    vp = a*VP/sqrt(gpp) - bp;


    c->prim[RHO] = rho;
    c->prim[PPP] = P;
    c->prim[URR] = vr;
    c->prim[UPP] = vp;
    c->prim[UZZ] = 0.0;

    for(i=sim_NUM_C(theSim); i<sim_NUM_Q(theSim); i++)
        c->prim[i] = q;
}

void cell_init_kick_calc(struct Cell *c, double r, double phi, double z,
                            struct Sim *theSim)
{
    int opt = sim_InitPar0(theSim);

    if(opt == 0)
        cell_init_kick_isothermal(c, r, phi, z, theSim);
    else
        printf("ERROR: cell_init_kick given bad option.\n");
}

void cell_single_init_kick(struct Cell *theCell, struct Sim *theSim,
                            int i, int j, int k)
{
    double rm = sim_FacePos(theSim,i-1,R_DIR);
    double rp = sim_FacePos(theSim,i,R_DIR);
    double r = 0.5*(rm+rp);
    double zm = sim_FacePos(theSim,k-1,Z_DIR);
    double zp = sim_FacePos(theSim,k,Z_DIR);
    double z = 0.5*(zm+zp);
    double t = theCell->tiph-.5*theCell->dphi;

    cell_init_kick_calc(theCell, r, t, z, theSim);
}

void cell_init_kick(struct Cell ***theCells,struct Sim *theSim,
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

                cell_init_kick_calc(&(theCells[k][i][j]), r, t, z, theSim);

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
