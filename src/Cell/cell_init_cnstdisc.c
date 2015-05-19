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
    double GAM = sim_GAMMALAW(theSim);
    double M = sim_GravM(theSim);
    double rho0 = sim_InitPar1(theSim);
    double T = sim_InitPar2(theSim);

    double rm = sim_FacePos(theSim,i-1,R_DIR);
    double rp = sim_FacePos(theSim,i,R_DIR);
    double r = 0.5*(rm+rp);

    double H = sqrt((1-3*M/r)*r*r*r*T/(M*(1.0+GAM*T/(GAM-1.0))));

    if(sim_Background(theSim) != GRDISC)
    {
        rho = rho0*H;
        Pp = rho0 * T *H;
    }
    else
    {
        rho = rho0;
        Pp = T;
    }
    
    vr = 0.0; 
    vp = sqrt(M/(r*r*r));

    theCell->prim[RHO] = rho;
    theCell->prim[PPP] = Pp;
    theCell->prim[URR] = vr;
    theCell->prim[UPP] = vp;
    theCell->prim[UZZ] = 0.0;
}

void cell_init_cnstdisc(struct Cell ***theCells,struct Sim *theSim,struct MPIsetup * theMPIsetup)
{

    double rho, Pp, vr, vp;
    double GAM = sim_GAMMALAW(theSim);
    double M = sim_GravM(theSim);
    double rho0 = sim_InitPar1(theSim);
    double T = sim_InitPar2(theSim);

    int i, j, k;
    for (k = 0; k < sim_N(theSim,Z_DIR); k++) 
    {
        for (i = 0; i < sim_N(theSim,R_DIR); i++) 
        {
            double rm = sim_FacePos(theSim,i-1,R_DIR);
            double rp = sim_FacePos(theSim,i,R_DIR);
            double r = 0.5*(rm+rp);
            
            double H = sqrt((1-3*M/r)*r*r*r*T/(M*(1.0+GAM*T/(GAM-1.0))));

            if(sim_Background(theSim) != GRDISC)
            {
                rho = rho0*H;
                Pp = rho0 * T *H;
            }
            else
            {
                rho = rho0;
                Pp = T;
            }
            
            vr = 0.0; 
            vp = sqrt(M/(r*r*r));

            for (j = 0; j < sim_N_p(theSim,i); j++) 
            {
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
