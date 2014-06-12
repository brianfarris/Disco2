#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Sim.h"
#include "../Headers/Face.h"
#include "../Headers/GravMass.h"
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
    double u0 = up/vp;
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
        double x;
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
    int test_num = sim_InitPar0(theSim);

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