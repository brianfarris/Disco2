#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Sim.h"
#include "../Headers/Face.h"
#include "../Headers/GravMass.h"
#include "../Headers/header.h"

// Constant Density Disc
void cell_single_init_cnstdisc(struct Cell *theCell, struct Sim *theSim,int i,int j,int k)
{
    double rho, Pp, vr, vp;
    double GAMMALAW = sim_GAMMALAW(theSim);
    double M = sim_GravM(theSim);
    double R0 = sim_InitPar1(theSim);
    double cs20 = sim_InitPar2(theSim);

    double rm = sim_FacePos(theSim,i-1,R_DIR);
    double rp = sim_FacePos(theSim,i,R_DIR);
    double r = 0.5*(rm+rp);
    double zm = sim_FacePos(theSim,k-1,Z_DIR);
    double zp = sim_FacePos(theSim,k,Z_DIR);
    double z = 0.5*(zm+zp);
    double t = theCell->tiph-.5*theCell->dphi;

    rho = R0;
    Pp = R0 * cs20;
    vr = 0.0; 
    vp = sqrt(M/(r*r*r));

    theCell->prim[RHO] = rho;
    theCell->prim[PPP] = Pp;
    theCell->prim[URR] = vr;
    theCell->prim[UPP] = vp;
    theCell->prim[UZZ] = 0.0;
    theCell->divB = 0.0;
    theCell->GradPsi[0] = 0.0;
    theCell->GradPsi[1] = 0.0;

}

void cell_init_cnstdisc(struct Cell ***theCells,struct Sim *theSim,struct MPIsetup * theMPIsetup)
{

    double rho, Pp, vr, vp;
    double GAMMALAW = sim_GAMMALAW(theSim);
    double M = sim_GravM(theSim);
    double R0 = sim_InitPar1(theSim);
    double cs20 = sim_InitPar2(theSim);

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
         
            rho = R0;
            Pp = R0 * cs20;
            vr = 0.0; 
            vp = sqrt(M/(r*r*r));

            for (j = 0; j < sim_N_p(theSim,i); j++) 
            {
                double t = theCells[k][i][j].tiph-.5*theCells[k][i][j].dphi;
             
                theCells[k][i][j].prim[RHO] = rho;
                theCells[k][i][j].prim[PPP] = Pp;
                theCells[k][i][j].prim[URR] = vr;
                theCells[k][i][j].prim[UPP] = vp;
                theCells[k][i][j].prim[UZZ] = 0.0;
                theCells[k][i][j].divB = 0.0;
                theCells[k][i][j].GradPsi[0] = 0.0;
                theCells[k][i][j].GradPsi[1] = 0.0;
                theCells[k][i][j].GradPsi[2] = 0.0;
            }
        }
    }
}
