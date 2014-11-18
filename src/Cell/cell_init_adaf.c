#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Sim.h"
#include "../Headers/Face.h"
#include "../Headers/GravMass.h"
#include "../Headers/header.h"

//ADAF: Advection-Dominated Accretion Flow

void cell_init_adaf_nocool(struct Cell *c, double r, double phi, double z, struct Sim *theSim)
{
    double GAM = sim_GAMMALAW(theSim);
    double M = sim_GravM(theSim);
    double alpha = sim_AlphaVisc(theSim);

    double rho, vr, vp, Pp;
    double rho0, vr0, vp0, P0, r0;

    rho0 = 1.0;
    r0 = sim_InitPar1(theSim);
    P0 = sim_InitPar2(theSim);
    vr0 = sim_InitPar3(theSim);
    vp0 = sim_InitPar4(theSim);

    rho = rho0 * pow(r/r0, -0.5);
    Pp = P0 * pow(r/r0, -1.5);
    if(r > 2*M)
    {
        vr = vr0 * pow((r+4*M)/(r0+4*M), -0.5);
        vp = vp0 * pow((r+4*M)/(r0+4*M), -1.5);
    }
    else
    {
        vr = vr0 * pow((6*M)/(r0+4*M), -0.5);
        vp = vp0 * pow((6*M)/(r0+4*M), -1.5);
    }


    if(sim_Metric(theSim) == SCHWARZSCHILD_KS)
    {
        double fac = 1.0/(1.0 + 2*M*vr/(r-2*M));
        vr *= fac;
        vp *= fac;

        if(1.0-2*M/r-4*M*vr/r-vr*vr*(1+2*M/r)-vp*vp*r*r < 0.0)
        {
            vp = 0.0;
            //vr = -2*M/(r+2*M);
            //if(r < 2*M)
            //    vr = (r/(2*M))*((r-2*M)/(r+2*M)) + (1-r/(2*M))*(-1.0);
            //else
                vr = 0.9*(r-2*M)/(r+2*M) + 0.1*(-1.0);
        }
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

void cell_init_adaf_calc(struct Cell *c, double r, double phi, double z, struct Sim *theSim)
{
    int disc_num = sim_InitPar0(theSim);

    if(disc_num == 0)
        cell_init_adaf_nocool(c, r, phi, z, theSim);
    else
        printf("ERROR: cell_init_accdisc given bad option.\n");
}

void cell_single_init_adaf(struct Cell *theCell, struct Sim *theSim,int i,int j,int k)
{
    double rm = sim_FacePos(theSim,i-1,R_DIR);
    double rp = sim_FacePos(theSim,i,R_DIR);
    double r = 0.5*(rm+rp);
    double zm = sim_FacePos(theSim,k-1,Z_DIR);
    double zp = sim_FacePos(theSim,k,Z_DIR);
    double z = 0.5*(zm+zp);
    double t = theCell->tiph-.5*theCell->dphi;

    cell_init_adaf_calc(theCell, r, t, z, theSim);
}

void cell_init_adaf(struct Cell ***theCells,struct Sim *theSim,struct MPIsetup * theMPIsetup)
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
             
                cell_init_adaf_calc(&(theCells[k][i][j]), r, t, z, theSim);
                
                if(PRINTTOOMUCH)
                {
                    printf("(%d,%d,%d) = (%.12lg, %.12lg, %.12lg): (%.12lg, %.12lg, %.12lg, %.12lg, %.12lg)\n", i,j,k,r,t,z,theCells[k][i][j].prim[RHO],theCells[k][i][j].prim[URR],theCells[k][i][j].prim[UPP],theCells[k][i][j].prim[UZZ],theCells[k][i][j].prim[PPP]);
                }
            }
        }
    }
}
