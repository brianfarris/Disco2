#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Sim.h"
#include "../Headers/Face.h"
#include "../Headers/GravMass.h"
#include "../Headers/header.h"

void cell_single_init_isentrope(struct Cell *theCell, struct Sim *theSim,int i,int j,int k){
    double rho_ref  = 1.0;
    double P_ref  = 100.0;
    double v_ref = 0.0;
    double GAMMALAW = sim_GAMMALAW(theSim);
    double L = 0.3;
    double a = 1.0;
    double x0 = 0.5;
    double rho, P, v;

    double rm = sim_FacePos(theSim,i-1,R_DIR);
    double rp = sim_FacePos(theSim,i,R_DIR);
    double r = 0.5*(rm+rp);
    double zm = sim_FacePos(theSim,k-1,Z_DIR);
    double zp = sim_FacePos(theSim,k,Z_DIR);
    double z = 0.5*(zm+zp);
    double t = theCell->tiph-.5*theCell->dphi;
    double x = (r*cos(t) - x0) / L;
    //x = (z-x0) / L;

    if(fabs(x) < 1.0)
    {
        double cs, cs_ref, Jm, Js;
        rho = rho_ref*(1.0 + a*(x*x-1.0)*(x*x-1.0)*(x*x-1.0)*(x*x-1.0));
        P = P_ref * pow(rho / rho_ref, GAMMALAW);
        
        cs = sqrt(GAMMALAW * P / (rho + GAMMALAW * P /(GAMMALAW-1.0)));
        cs_ref = sqrt(GAMMALAW * P_ref / (rho_ref + GAMMALAW * P_ref /(GAMMALAW-1.0)));
        Jm = 0.5 * log((1.0+v_ref)/(1.0-v_ref)) - log((sqrt(GAMMALAW-1)+cs_ref)/(sqrt(GAMMALAW-1)-cs_ref)) / sqrt(GAMMALAW-1);
        Js = exp(2*Jm + 2*log((sqrt(GAMMALAW-1)+cs)/(sqrt(GAMMALAW-1)-cs))/sqrt(GAMMALAW-1));
        
        v = (Js - 1.0) / (Js + 1.0);
    }
    else
    {
        rho = rho_ref;
        P = P_ref;
        v = v_ref;
    }

    theCell->prim[RHO] = rho;
    theCell->prim[PPP] = P;
    theCell->prim[URR] = v * cos(t);
    theCell->prim[UPP] = -v * sin(t) / r;
    theCell->prim[UZZ] = 0.0;
    theCell->divB = 0.0;
    theCell->GradPsi[0] = 0.0;
    theCell->GradPsi[1] = 0.0;

  //TODO: Not sure what this is for.  Ask someone if important.
  //if(sim_NUM_C(theSim)<sim_NUM_Q(theSim)) theCell->prim[sim_NUM_C(theSim)] = Qq;
}

void cell_init_isentrope(struct Cell ***theCells,struct Sim *theSim,struct MPIsetup * theMPIsetup)
{
    double rho_ref  = 1.0;
    double P_ref  = 100.0;
    double v_ref  = 0.0;
    double GAMMALAW = sim_GAMMALAW(theSim);
    double L = 0.3;
    double a = 1.0;
    double x0 = 0.5;
    double K = P_ref / pow(rho_ref, GAMMALAW);
    double rho, P, v;

    int i, j, k;
    for (k = 0; k < sim_N(theSim,Z_DIR); k++) 
    {
        for (i = 0; i < sim_N(theSim,R_DIR); i++) 
        {
            double rm = sim_FacePos(theSim,i-1,R_DIR);
            double rp = sim_FacePos(theSim,i,R_DIR);
            double r = 0.5*(rm+rp);
            double zm = sim_FacePos(theSim,k-1,Z_DIR);
            double zp = sim_FacePos(theSim,k,Z_DIR);
            double z = 0.5*(zm+zp);
            for (j = 0; j < sim_N_p(theSim,i); j++) 
            {
                double t = theCells[k][i][j].tiph-.5*theCells[k][i][j].dphi;
                double x = (r*cos(t) - x0) / L;
                //x = (z-x0) / L;

                if(fabs(x) < 1.0)
                {
                    double cs, cs_ref, Jm, Js;
                    rho = rho_ref*(1.0 + a*(x*x-1.0)*(x*x-1.0)*(x*x-1.0)*(x*x-1.0));
                    P = P_ref * pow(rho / rho_ref, GAMMALAW);
                    
                    cs = sqrt(GAMMALAW * P / (rho + GAMMALAW * P /(GAMMALAW-1.0)));
                    cs_ref = sqrt(GAMMALAW * P_ref / (rho_ref + GAMMALAW * P_ref /(GAMMALAW-1.0)));
                    Jm = 0.5 * log((1.0+v_ref)/(1.0-v_ref)) - log((sqrt(GAMMALAW-1)+cs_ref)/(sqrt(GAMMALAW-1)-cs_ref)) / sqrt(GAMMALAW-1);
                    Js = exp(2*Jm + 2*log((sqrt(GAMMALAW-1)+cs)/(sqrt(GAMMALAW-1)-cs))/sqrt(GAMMALAW-1));
                    
                    v = (Js - 1.0) / (Js + 1.0);
                }
                else
                {
                    rho = rho_ref;
                    P = P_ref;
                    v = v_ref;
                }

                theCells[k][i][j].prim[RHO] = rho;
                theCells[k][i][j].prim[PPP] = P;
                theCells[k][i][j].prim[URR] = v * cos(t);
                theCells[k][i][j].prim[UPP] = -v * sin(t) / r;
                theCells[k][i][j].prim[UZZ] = 0.0;
                theCells[k][i][j].divB = 0.0;
                theCells[k][i][j].GradPsi[0] = 0.0;
                theCells[k][i][j].GradPsi[1] = 0.0;
                theCells[k][i][j].GradPsi[2] = 0.0;

            }
        }
    }
}
