#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Sim.h"
#include "../Headers/Face.h"
#include "../Headers/GravMass.h"
#include "../Headers/header.h"

//Accreting Disc
//
//This is meant to be a generic initial condition for relativistic discs.  It initializes with a
//constant density and pressure, keplerian orbital velocity, and constant accretion-rate radial
//velocity.

void cell_init_accdisc_flat(struct Cell *c, double r, double phi, double z, struct Sim *theSim)
{
    double Mdot = sim_InitPar1(theSim);
    double rho = sim_InitPar2(theSim);
    double Pp = sim_InitPar3(theSim);
    double GAM = sim_GAMMALAW(theSim);
    double M = sim_GravM(theSim);

    double vr;
    double vp;
        
    double ur = -Mdot/(2*M_PI*rho*r);
    //vp = exp(-1.0/(r/M-2)) * sqrt(M/(r*r*r));
    vp = exp(-0.5/(r/M-3)) * sqrt(M/(r*r*r));

    if(sim_Metric(theSim) == SCHWARZSCHILD_SC)
    {
        //Only valid outside r = 2*M
        double u0 = sqrt((1+ur*ur/(1-2*M/r))/(1-2*M/r-r*r*vp*vp));
        vr = ur/u0;
    }
    else if(sim_Metric(theSim) == SCHWARZSCHILD_KS)
    {
        if(r <= 3*M)
            vp = 0.0;
        if(r <=2*M && ur*ur +1.0 < 2*M/r)
            ur = -1.01 * sqrt(2*M/r - 1.0);
        double u0 = (2*M*ur/r + sqrt((1-(1+2*M/r)*r*r*vp*vp)*ur*ur + 1.0-2*M/r-r*r*vp*vp)) / (1-2*M/r-r*r*vp*vp);
        vr = ur/u0;
    }
    
    c->prim[RHO] = rho;
    c->prim[PPP] = Pp;
    c->prim[URR] = vr;
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

void cell_init_accdisc_calc(struct Cell *c, double r, double phi, double z, struct Sim *theSim)
{
    int disc_num = sim_InitPar0(theSim);

    if(disc_num == 0)
        cell_init_accdisc_flat(c, r, phi, z, theSim);
    else
        printf("ERROR: cell_init_accdisc given bad option.\n");
}

void cell_single_init_accdisc(struct Cell *theCell, struct Sim *theSim,int i,int j,int k)
{
    double rm = sim_FacePos(theSim,i-1,R_DIR);
    double rp = sim_FacePos(theSim,i,R_DIR);
    double r = 0.5*(rm+rp);
    double zm = sim_FacePos(theSim,k-1,Z_DIR);
    double zp = sim_FacePos(theSim,k,Z_DIR);
    double z = 0.5*(zm+zp);
    double t = theCell->tiph-.5*theCell->dphi;

    cell_init_accdisc_calc(theCell, r, t, z, theSim);
}

void cell_init_accdisc(struct Cell ***theCells,struct Sim *theSim,struct MPIsetup * theMPIsetup)
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
             
                cell_init_accdisc_calc(&(theCells[k][i][j]), r, t, z, theSim);
                
                if(PRINTTOOMUCH)
                {
                    printf("(%d,%d,%d) = (%.12lg, %.12lg, %.12lg): (%.12lg, %.12lg, %.12lg, %.12lg, %.12lg)\n", i,j,k,r,t,z,theCells[k][i][j].prim[RHO],theCells[k][i][j].prim[URR],theCells[k][i][j].prim[UPP],theCells[k][i][j].prim[UZZ],theCells[k][i][j].prim[PPP]);
                }
            }
        }
    }
}
