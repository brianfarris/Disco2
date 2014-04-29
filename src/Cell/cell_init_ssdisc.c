#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Sim.h"
#include "../Headers/Face.h"
#include "../Headers/GravMass.h"
#include "../Headers/header.h"

// Shakura-Sunyaev Disc
void cell_single_init_ssdisc(struct Cell *theCell, struct Sim *theSim,int i,int j,int k)
{
    double rho, Pp, vr, vp;
    double GAMMALAW = sim_GAMMALAW(theSim);
    double M = sim_GravM(theSim);
    double Rs = sim_InitPar1(theSim);
    double Mdot = sim_InitPar2(theSim);
    double alpha = sim_AlphaVisc(theSim);
    double q0 = sim_CoolFac(theSim);

    double rm = sim_FacePos(theSim,i-1,R_DIR);
    double rp = sim_FacePos(theSim,i,R_DIR);
    double r = 0.5*(rm+rp);
    double zm = sim_FacePos(theSim,k-1,Z_DIR);
    double zp = sim_FacePos(theSim,k,Z_DIR);
    double z = 0.5*(zm+zp);
    double t = theCell->tiph-.5*theCell->dphi;

    double f = 1.0 - sqrt(Rs/r);

    double rho0 = pow(2.0,0.75) * pow(M_PI,-0.55) / 3.0;
    double vr0 = pow(2.0,-1.5) * pow(M_PI,-0.3) * 3;
    double P0 = pow(2.0,0.25) * pow(M_PI,-0.85) / 3.0;

    rho = rho0 * pow(alpha,-0.7) * pow(Mdot,0.55) * pow(M,0.625) * pow(r,-1.875) * pow(f,0.55) * pow(0.75*q0,0.15);
    vr = -vr0 * pow(alpha,0.8) * pow(Mdot,0.3) * pow(M,-0.25) * pow(r,-0.25) * pow(f,-0.7) * pow(0.75*q0,-0.1);
    vp = sqrt(M/(r*r*r));
    Pp = P0 * pow(alpha,-0.9) * pow(Mdot,0.85) * pow(M,0.875) * pow(r,-2.625) * pow(f,0.85) * pow(0.75*q0,0.05);

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

void cell_init_ssdisc(struct Cell ***theCells,struct Sim *theSim,struct MPIsetup * theMPIsetup)
{

    double rho, Pp, vr, vp;
    double GAMMALAW = sim_GAMMALAW(theSim);
    double M = sim_GravM(theSim);
    double Rs = sim_InitPar1(theSim);
    double Mdot = sim_InitPar2(theSim);
    double alpha = sim_AlphaVisc(theSim);
    double q0 = sim_CoolFac(theSim);

    double rho0 = pow(2.0,0.75) * pow(M_PI,-0.55) / 3.0;
    double vr0 = pow(2.0,-1.5) * pow(M_PI,-0.3) * 3;
    double P0 = pow(2.0,0.25) * pow(M_PI,-0.85) / 3.0;
    double f;

    rho0 *= pow(alpha,-0.7) * pow(Mdot,0.55) * pow(M,0.625) * pow(0.75*q0,0.15);
    vr0 *= pow(alpha,0.8) * pow(Mdot,0.3) * pow(M,-0.25) * pow(0.75*q0,-0.1);
    P0 *= pow(alpha,-0.9) * pow(Mdot,0.85) * pow(M,0.875) * pow(0.75*q0,0.05);

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
            f = 1.0 - sqrt(Rs/r);

            rho = rho0 * pow(r,-1.875) * pow(f,0.55);
            vr = -vr0 * pow(r,-0.25) * pow(f,-0.7);
            vp = sqrt(M/(r*r*r));
            Pp = P0 * pow(r,-2.625) * pow(f,0.85);

            if(sim_Metric(theSim) == SCHWARZSCHILD_KS)
            {
                vr  = vr * (1.0-M/r) / (1.0-M/r+M/r*vr);
                if(vr > (1.0-M/r) / (1.0+M/r))
                    vr = 0.9 * (1.0-M/r) / (1.0+M/r) - 0.1*vr;
            }

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
