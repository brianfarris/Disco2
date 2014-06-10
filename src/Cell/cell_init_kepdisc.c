#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Sim.h"
#include "../Headers/Face.h"
#include "../Headers/GravMass.h"
#include "../Headers/header.h"

// Keplerian Disc
void cell_single_init_kepdisc(struct Cell *theCell, struct Sim *theSim,int i,int j,int k)
{
    double rho, Pp, vr, vp;
    double GAMMALAW = sim_GAMMALAW(theSim);
    double B = sim_InitPar1(theSim);
    double mach0 = sim_InitPar2(theSim);
    double r0 = sim_InitPar3(theSim);
    double dr = sim_InitPar4(theSim);
    double nu = fabs(sim_AlphaVisc(theSim));
    double M = sim_GravM(theSim);

    double rm = sim_FacePos(theSim,i-1,R_DIR);
    double rp = sim_FacePos(theSim,i,R_DIR);
    double r = 0.5*(rm+rp);
    double zm = sim_FacePos(theSim,k-1,Z_DIR);
    double zp = sim_FacePos(theSim,k,Z_DIR);
    double z = 0.5*(zm+zp);
    double t = theCell->tiph-.5*theCell->dphi;

    rho = 1.0 + B/sqrt(r) + exp(-(r-r0)*(r-r0)/(2*dr*dr));
    Pp = (1.0 + B/sqrt(r0)) * M / (r0 * GAMMALAW * mach0*mach0);
    vr = - 3*nu*(1.0+(1.0-r*(r-r0)/(dr*dr))*exp(-(r-r0)*(r-r0)/(2*dr*dr))) / (2*r*rho);
    vp = sqrt(M/(r*r*r));
    if(dr < 0)
    {
        rho = 1.0 + B/sqrt(r);
        vr = -3*nu/(2*r*rho);
    }
   
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

void cell_init_kepdisc(struct Cell ***theCells,struct Sim *theSim,struct MPIsetup * theMPIsetup)
{

    double rho, Pp, vr, vp;
    double GAMMALAW = sim_GAMMALAW(theSim);
    double B = sim_InitPar1(theSim);
    double mach0 = sim_InitPar2(theSim);
    double r0 = sim_InitPar3(theSim);
    double dr = sim_InitPar4(theSim);
    double M = sim_GravM(theSim);
    double nu = fabs(sim_AlphaVisc(theSim));

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

            rho = 1.0 + B/sqrt(r) + exp(-(r-r0)*(r-r0)/(2*dr*dr));
            Pp = (1.0 + B/sqrt(r0)) * M / (r0 * GAMMALAW * mach0*mach0);
            vr = - 3*nu*(1.0+(1.0-r*(r-r0)/(dr*dr))*exp(-(r-r0)*(r-r0)/(2*dr*dr))) / (2*r*rho);
            vp = sqrt(M/(r*r*r));
            if(dr < 0)
            {
                rho = 1.0 + B/sqrt(r);
                vr = -3*nu/(2*r*rho);
            }

            for (j = 0; j < sim_N_p(theSim,i); j++) 
            {
                double t = theCells[k][i][j].tiph-.5*theCells[k][i][j].dphi;
            
                if(sim_Metric(theSim) == SCHWARZSCHILD_KS)
                {
                    vr  = vr * (1.0-M/r) / (1.0-M/r+M/r*vr);
                    if(vr > (1.0-M/r) / (1.0+M/r))
                        vr = 0.9 * (1.0-M/r) / (1.0+M/r) - 0.1*vr;
                }
             
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
                    printf("(%d,%d,%d) = (%.12lg, %.12lg, %.12lg): (%.12lg, %.12lg, %.12lg, %.12lg, %.12lg)\n", i,j,k,r,t,z,rho,vr,vp,theCells[k][i][j].prim[UZZ],Pp);
                }
            }
        }
    }
}
