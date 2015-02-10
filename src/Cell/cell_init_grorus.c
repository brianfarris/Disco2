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

//GRorus: GR(T)orus solution.

void cell_init_grorus_newt(struct Cell *c, double r, double phi, double z,
                            struct Sim *theSim)
{
    int i;
    double M = sim_GravM(theSim);
    double GAM = sim_GAMMALAW(theSim);

    double rho, P, vr, vp, vz;
    double C, j, cs20, cs2, rho_floor, cs2_floor;

    C = sim_InitPar1(theSim);
    j = sim_InitPar2(theSim);
    cs20 = sim_InitPar3(theSim);

    cs2 = (GAM-1.0)/GAM * (C + M/sqrt(r*r+z*z) - 0.5*j*j/(r*r));

    rho_floor = sim_RHO_FLOOR(theSim);
    cs2_floor = cs20 * pow(rho_floor,GAM-1.0);

    if(cs2 > cs2_floor)
        rho = pow(cs2/cs20, 1.0/(GAM-1.0));
    else
    {
        rho = rho_floor;
        cs2 = cs2_floor;
    }

    P = cs2*rho;
    vr = 0.0;
    vp = j/(r*r);
    vz = 0.0;

    c->prim[RHO] = rho;
    c->prim[PPP] = P;
    c->prim[URR] = vr;
    c->prim[UPP] = vp;
    c->prim[UZZ] = vz;

    for(i=sim_NUM_C(theSim); i<sim_NUM_Q(theSim); i++)
    {
        if(r*cos(phi) < 0)
            c->prim[i] = 0.0;
        else
            c->prim[i] = 1.0;
    }
}

void cell_init_grorus_calc(struct Cell *c, double r, double phi, double z,
                            struct Sim *theSim)
{
    int choice = sim_InitPar0(theSim);

    if(choice == 0)
        cell_init_grorus_newt(c, r, phi, z, theSim);
    else
        printf("ERROR: cell_init_grorus given bad option.\n");
}

void cell_single_init_grorus(struct Cell *theCell, struct Sim *theSim,
                            int i, int j, int k)
{
    double rm = sim_FacePos(theSim,i-1,R_DIR);
    double rp = sim_FacePos(theSim,i,R_DIR);
    double r = 0.5*(rm+rp);
    double zm = sim_FacePos(theSim,k-1,Z_DIR);
    double zp = sim_FacePos(theSim,k,Z_DIR);
    double z = 0.5*(zm+zp);
    double t = theCell->tiph-.5*theCell->dphi;

    cell_init_grorus_calc(theCell, r, t, z, theSim);
}

void cell_init_grorus(struct Cell ***theCells,struct Sim *theSim,
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

                cell_init_grorus_calc(&(theCells[k][i][j]), r, t, z, theSim);

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
