#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Sim.h"
#include "../Headers/Face.h"
#include "../Headers/GravMass.h"
#include "../Headers/header.h"

// Cartesian Shear Test
void cell_single_init_cartshear(struct Cell *theCell, struct Sim *theSim,int i,int j,int k)
{
    double rho, Pp, vr, vp;
    double GAMMALAW = sim_GAMMALAW(theSim);
    double v0 = sim_InitPar1(theSim);
    double cs20 = sim_InitPar2(theSim);
    double t0 = sim_InitPar3(theSim);
    double x0 = sim_InitPar4(theSim);
    double alpha = fabs(sim_AlphaVisc(theSim));

    double rm = sim_FacePos(theSim,i-1,R_DIR);
    double rp = sim_FacePos(theSim,i,R_DIR);
    double r = 0.5*(rm+rp);
    double zm = sim_FacePos(theSim,k-1,Z_DIR);
    double zp = sim_FacePos(theSim,k,Z_DIR);
    double z = 0.5*(zm+zp);
    double t = theCell->tiph-.5*theCell->dphi;

    double rho0 = 1.0;
    double P0 = cs20 * rho0;
    double dx = r * cos(t) - x0;
    double vy = v0 * exp(-dx*dx/(4.0*alpha*t0)) / sqrt(4.0*M_PI*alpha*t0);

    rho = rho0;
    Pp = P0;
    vr = vy * sin(t); 
    vp = vy * cos(t) / r;

    theCell->prim[RHO] = rho;
    theCell->prim[PPP] = Pp;
    theCell->prim[URR] = vr;
    theCell->prim[UPP] = vp;
    theCell->prim[UZZ] = 0.0;
    theCell->divB = 0.0;
    theCell->GradPsi[0] = 0.0;
    theCell->GradPsi[1] = 0.0;

  //TODO: Not sure what this is for.  Ask someone if important.
  //if(sim_NUM_C(theSim)<sim_NUM_Q(theSim)) theCell->prim[sim_NUM_C(theSim)] = Qq;
}

void cell_init_cartshear(struct Cell ***theCells,struct Sim *theSim,struct MPIsetup * theMPIsetup)
{

    double rho, Pp, vr, vp;
    double GAMMALAW = sim_GAMMALAW(theSim);
    double M = sim_GravM(theSim);
    double v0 = sim_InitPar1(theSim);
    double cs20 = sim_InitPar2(theSim);
    double t0 = sim_InitPar3(theSim);
    double x0 = sim_InitPar4(theSim);
    double alpha = fabs(sim_AlphaVisc(theSim));

    double rho0 = 1.0;
    double P0 = cs20 * rho0;
    
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

            rho = rho0;
            Pp = P0;

            if(sim_Metric(theSim) == SCHWARZSCHILD_KS)
            {
                vr  = vr * (1.0-M/r) / (1.0-M/r+M/r*vr);
                if(vr > (1.0-M/r) / (1.0+M/r))
                    vr = 0.9 * (1.0-M/r) / (1.0+M/r) - 0.1*vr;
            }

            for (j = 0; j < sim_N_p(theSim,i); j++) 
            {
                double t = theCells[k][i][j].tiph-.5*theCells[k][i][j].dphi;
                double dx = r*cos(t) - x0;
                double vy = v0 * exp(-dx*dx/(4.0*alpha*t0)) / sqrt(4.0*M_PI*alpha*t0);
                vr = vy * sin(t); 
                vp = vy * cos(t) / r;
             
                theCells[k][i][j].prim[RHO] = rho;
                theCells[k][i][j].prim[PPP] = Pp;
                theCells[k][i][j].prim[URR] = vr;
                theCells[k][i][j].prim[UPP] = vp;
                theCells[k][i][j].prim[UZZ] = 0.0;
                theCells[k][i][j].divB = 0.0;
                theCells[k][i][j].GradPsi[0] = 0.0;
                theCells[k][i][j].GradPsi[1] = 0.0;
                theCells[k][i][j].GradPsi[2] = 0.0;

                if(PRINTTOOMUCH)
                {
                    printf("(%d,%d,%d) = (%.12lg, %.12lg, %.12lg): (%.12lg, %.12lg, %.12lg, %.12lg)\n", i,j,k,r,t,z,rho,vr,vp,Pp);
                }
            }
        }
    }
}
