#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Sim.h"
#include "../Headers/Face.h"
#include "../Headers/GravMass.h"
#include "../Headers/header.h"

//Pressure supported hydrostatic equilibrium for the Schwarzschild Metric.

void cell_single_init_equil1(struct Cell *theCell, struct Sim *theSim,int i,int j,int k)
{
    double rho, Pp;
    double GAMMALAW = sim_GAMMALAW(theSim);
    double RG = sim_GravRadius(theSim);
    double rho0 = sim_InitPar1(theSim);
    double P0 = sim_InitPar2(theSim);
    
    double rm = sim_FacePos(theSim,i-1,R_DIR);
    double rp = sim_FacePos(theSim,i,R_DIR);
    double r = 0.5*(rm+rp);
    double zm = sim_FacePos(theSim,k-1,Z_DIR);
    double zp = sim_FacePos(theSim,k,Z_DIR);
    double z = 0.5*(zm+zp);
    double t = theCell->tiph-.5*theCell->dphi;

    double a = (GAMMALAW-1.0)/GAMMALAW;
    
    rho = rho0;
    Pp = (P0+rho0*a)*pow(1.0-RG/r, -0.5/a) - rho0*a;

    theCell->prim[RHO] = rho;
    theCell->prim[PPP] = Pp;
    theCell->prim[URR] = 0.0;
    theCell->prim[UPP] = 0.0;
    theCell->prim[UZZ] = 0.0;
    theCell->divB = 0.0;
    theCell->GradPsi[0] = 0.0;
    theCell->GradPsi[1] = 0.0;

  //TODO: Not sure what this is for.  Ask someone if important.
  //if(sim_NUM_C(theSim)<sim_NUM_Q(theSim)) theCell->prim[sim_NUM_C(theSim)] = Qq;
}

void cell_init_equil1(struct Cell ***theCells,struct Sim *theSim,struct MPIsetup * theMPIsetup)
{

    double rho, Pp;
    double GAMMALAW = sim_GAMMALAW(theSim);
    double RG = sim_GravRadius(theSim);
    double rho0 = sim_InitPar1(theSim);
    double P0 = sim_InitPar2(theSim);

    double a = (GAMMALAW-1.0)/GAMMALAW;
    rho = rho0;

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
            
            Pp = (P0+rho0*a)*pow(1.0-RG/r, -0.5/a) - rho0*a;

            for (j = 0; j < sim_N_p(theSim,i); j++) 
            {
                double t = theCells[k][i][j].tiph-.5*theCells[k][i][j].dphi;
             
                theCells[k][i][j].prim[RHO] = rho;
                theCells[k][i][j].prim[PPP] = Pp;
                theCells[k][i][j].prim[URR] = 0.0;
                theCells[k][i][j].prim[UPP] = 0.0;
                theCells[k][i][j].prim[UZZ] = 0.0;
                theCells[k][i][j].divB = 0.0;
                theCells[k][i][j].GradPsi[0] = 0.0;
                theCells[k][i][j].GradPsi[1] = 0.0;
                theCells[k][i][j].GradPsi[2] = 0.0;
            }
        }
    }
}
