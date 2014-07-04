#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Sim.h"
#include "../Headers/Face.h"
#include "../Headers/GravMass.h"
#include "../Headers/header.h"

//Shakura-Sunyaev Disc

//Thompson Cooling
void cell_init_ssdisc_thompson(struct Cell *c, double r, double phi, double z, struct Sim *theSim)
{
    double mach0 = sim_InitPar1(theSim);
    double r0 = sim_InitPar2(theSim);
    double DR = sim_InitPar3(theSim);
    double GAM = sim_GAMMALAW(theSim);
    double M = sim_GravM(theSim);
    double alpha = sim_AlphaVisc(theSim);
    double q0 = sim_CoolFac(theSim);

    if(q0 == 0.0)
        q0 = 1.0;

    double rho0 = sqrt(2*q0/alpha)*pow(M*M*M*M*M/(GAM*GAM*GAM*GAM*GAM*GAM*GAM*r0*r0*r0), 0.25) / (3.0*mach0*mach0*mach0);
    double P0 = M*rho0/(GAM*r0*mach0*mach0);
    double vr0 = -3*alpha * sqrt(M/(GAM*r0)) / (mach0*mach0);
    
    double rho = rho0 * pow(r0/r, 0.6);
    double P = P0 * pow(r0/r, 1.5);
    double vr = vr0 * pow(r0/r, 0.4);
    double vp = sqrt(M/(r*r*r));

    if(DR > 0.0 && r < r0)
    {
        double prof = exp(-(r-r0)*(r-r0)/(2*DR*DR));
        rho *= prof;
        P *= prof;
    }

    //TODO: Make General
    //if(r*vp >= sqrt(1-2*M/r-vr*vr/(1-2*M/r)))
    if(r <= 6*M)
        vp = 0.5*sqrt(1-2*M/r-vr*vr/(1-2*M/r))/r;
    
    c->prim[RHO] = rho;
    c->prim[PPP] = P;
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

void cell_init_ssdisc_isotherm(struct Cell *c, double r, double phi, double z, struct Sim *theSim)
{
    double mach0 = sim_InitPar1(theSim);
    double r0 = sim_InitPar2(theSim);
    double DR = sim_InitPar3(theSim);
    double rho0 = sim_InitPar4(theSim);
    double GAM = sim_GAMMALAW(theSim);
    double M = sim_GravM(theSim);
    double alpha = sim_AlphaVisc(theSim);
    double q0 = sim_CoolFac(theSim);

    double P0  = rho0 * M / (sqrt(GAM) * r0 * mach0*mach0);
    double vr0 = - 3 * alpha * sqrt(M/(GAM*r0)) / (mach0*mach0);
    
    double rho = rho0 * pow(r/r0, -1.5);
    double P = P0 * pow(r/r0, -1.5);
    double vr = vr0 * sqrt(r/r0);
    double vp = sqrt(M/(r*r*r));

    if(DR > 0.0 && r < r0)
    {
        double prof = exp(-(r-r0)*(r-r0)/(2*DR*DR));
        rho *= prof;
        P *= prof;
    }
    
    c->prim[RHO] = rho;
    c->prim[PPP] = P;
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
void cell_init_ssdisc_calc(struct Cell *c, double r, double phi, double z, struct Sim *theSim)
{
    int disc_num = sim_InitPar0(theSim);

    if(disc_num == 0)
        cell_init_ssdisc_thompson(c, r, phi, z, theSim);
    else if(disc_num == 1)
        cell_init_ssdisc_isotherm(c, r, phi, z, theSim);
    else
        printf("ERROR: cell_init_ssdisc given bad option.\n");
}

void cell_single_init_ssdisc(struct Cell *theCell, struct Sim *theSim,int i,int j,int k)
{
    double rm = sim_FacePos(theSim,i-1,R_DIR);
    double rp = sim_FacePos(theSim,i,R_DIR);
    double r = 0.5*(rm+rp);
    double zm = sim_FacePos(theSim,k-1,Z_DIR);
    double zp = sim_FacePos(theSim,k,Z_DIR);
    double z = 0.5*(zm+zp);
    double t = theCell->tiph-.5*theCell->dphi;

    cell_init_ssdisc_calc(theCell, r, t, z, theSim);
}

void cell_init_ssdisc(struct Cell ***theCells,struct Sim *theSim,struct MPIsetup * theMPIsetup)
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
             
                cell_init_ssdisc_calc(&(theCells[k][i][j]), r, t, z, theSim);
                
                if(PRINTTOOMUCH)
                {
                    printf("(%d,%d,%d) = (%.12lg, %.12lg, %.12lg): (%.12lg, %.12lg, %.12lg, %.12lg, %.12lg)\n", i,j,k,r,t,z,theCells[k][i][j].prim[RHO],theCells[k][i][j].prim[URR],theCells[k][i][j].prim[UPP],theCells[k][i][j].prim[UZZ],theCells[k][i][j].prim[PPP]);
                }
            }
        }
    }
}
