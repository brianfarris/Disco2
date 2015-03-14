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

//Orbit: Orbitting a binary.

void cell_init_orbit_passive(struct Cell *c, double r, double phi, double z,
                            struct Sim *theSim)
{
    double rho, P, vr, vp, vz, q;
    double rho0, P0, R0;
    int i;

    rho0 = sim_InitPar1(theSim);
    P0 = sim_InitPar2(theSim);
    R0 = sim_InitPar3(theSim);

    rho = rho0;
    P = P0;
    vr = 0.0;
    vp = 0.0;
    vz = 0.0;

    if(r < R0)
        q = 1.0;
    else
        q = 0.0;
    //q = r<R0 ? 1.0 : 0.0;

    c->prim[RHO] = rho;
    c->prim[PPP] = P;
    c->prim[URR] = vr;
    c->prim[UPP] = vp;
    c->prim[UZZ] = vz;

    for(i=sim_NUM_C(theSim); i<sim_NUM_Q(theSim); i++)
        c->prim[i] = q;
}

void cell_init_orbit_calc(struct Cell *c, double r, double phi, double z,
                            struct Sim *theSim)
{
    int opt = sim_InitPar0(theSim);

    if(opt == 0)
        cell_init_orbit_passive(c, r, phi, z, theSim);
    else
        printf("ERROR: cell_init_orbit given bad option.\n");
}

void cell_single_init_orbit(struct Cell *theCell, struct Sim *theSim,
                            int i, int j, int k)
{
    double rm = sim_FacePos(theSim,i-1,R_DIR);
    double rp = sim_FacePos(theSim,i,R_DIR);
    double r = 0.5*(rm+rp);
    double zm = sim_FacePos(theSim,k-1,Z_DIR);
    double zp = sim_FacePos(theSim,k,Z_DIR);
    double z = 0.5*(zm+zp);
    double t = theCell->tiph-.5*theCell->dphi;

    cell_init_orbit_calc(theCell, r, t, z, theSim);
}

void cell_init_orbit(struct Cell ***theCells,struct Sim *theSim,
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

                cell_init_orbit_calc(&(theCells[k][i][j]), r, t, z, theSim);

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
