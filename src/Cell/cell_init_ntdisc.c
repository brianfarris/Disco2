#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Sim.h"
#include "../Headers/Face.h"
#include "../Headers/GravMass.h"
#include "../Headers/header.h"

//Novikov-Thorne Disc

//Thompson Cooling
void cell_init_ntdisc_thompson(struct Cell *c, double r, double phi, double z, struct Sim *theSim)
{
    double Mdot = sim_InitPar1(theSim);
    double r0 = sim_InitPar2(theSim);
    double DR = sim_InitPar3(theSim);
    double GAM = sim_GAMMALAW(theSim);
    double M = sim_GravM(theSim);
    double alpha = sim_AlphaVisc(theSim);
    double q0 = sim_CoolPar1(theSim); //TODO: Assumes CoolingType == COOL_BB_ES
    double rs = 6*M;

    if(q0 == 0.0)
        q0 = 1.0;

    double C = 1.0-3*M/r;
    double D = 1.0-2*M/r;
    double P = 1.0 - sqrt(rs/r) + sqrt(3*M/r)*(atanh(sqrt(3*M/r)) - atanh(sqrt(3*M/rs)));
    
    double rho0 = 1.0/3.0*pow(Mdot*Mdot*Mdot*q0/(4*M_PI*M_PI*M_PI*alpha*alpha*alpha*alpha*GAM*GAM), 0.2);
    double P0 = Mdot/(6*M_PI*alpha*sqrt(GAM));
    double vr0 = -3.0*pow(alpha*alpha*alpha*alpha*GAM*GAM*Mdot*Mdot/(8*M_PI*M_PI*q0), 0.2);

    double rho = rho0 * pow(M/(r*r*r), 0.2) * pow(C,0.6) * pow(D,-1.6) * pow(P,0.6);
    double Pp = P0 * sqrt(M/(r*r*r)) * sqrt(C) * P / (D*D); 
    double vr = vr0 * pow(M*r*r, -0.2) * pow(C,-0.1) * pow(D,1.6) * pow(P,-0.6);
    double vp = sqrt(M/(r*r*r));
    double vz = 0.0;

    if(r <= 1.1*rs || r <= 3*M)
    {
        double rr = 1.1*rs;
        C = 1.0-3*M/rr;
        D = 1.0-2*M/rr;
        P = 1.0 - sqrt(rs/rr) + sqrt(3*M/rr)*(atanh(sqrt(3*M/rr)) - atanh(sqrt(3*M/rs)));
        
        rho = (1.0-0.9*(rr-r)/rr) * rho0 * pow(M/(rr*rr*rr), 0.2) * pow(C,0.6) * pow(D,-1.6) * pow(P,0.6);
        Pp = (1.0-0.9*(rr-r)/rr) * P0 * sqrt(M/(rr*rr*rr)) * sqrt(C) * P / (D*D); 
        vr = vr0 * pow(M*rr*rr, -0.2) * pow(C,-0.1) * pow(D,1.6) * pow(P,-0.6);
        vp = 0.5*sqrt(1-2*M/r-vr*vr/(1-2*M/r))/r;
    }

    if(DR > 0.0 && r < r0)
    {
        double prof = exp(-(r-r0)*(r-r0)/(2*DR*DR));
        rho *= prof;
        Pp *= prof;
    }

    if(sim_Metric(theSim) == SCHWARZSCHILD_KS)
    {
        double R = sqrt(r*r+z*z);
        double timefac = 1.0 / (1.0 + (r*vr+z*vz)/(R*R/(2*M)-1.0));
        vr *= timefac;
        vp *= timefac;
        vz *= timefac;

        if(r <= rs || r <= 3*M)
        {
            double rr = 1.1*rs;
            C = 1.0-3*M/rr;
            D = 1.0-2*M/rr;
            P = 1.0 - sqrt(rs/rr) + sqrt(3*M/rr)*(atanh(sqrt(3*M/rr)) - atanh(sqrt(3*M/rs)));
            rho = 0.5 * rho0 * pow(M/(rr*rr*rr), 0.2) * pow(C,0.6) * pow(D,-1.6) * pow(P,0.6);
            Pp = 0.5 * P0 * sqrt(M/(rr*rr*rr)) * sqrt(C) * P / (D*D); 
            if (r < 3*M)
            {
                vr = -(2.0*M/r) / (1.0+2.0*M/r);
                vp = 0.0;
            }
            else
            {
                vr = -(rs-r) / (rs-3*M) * (2.0*M/r) / (1.0+2.0*M/r);
                vp = (r-3*M)/(rs-3*M) * sqrt(M/(rr*rr*rr));
            }

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

void cell_init_ntdisc_calc(struct Cell *c, double r, double phi, double z, struct Sim *theSim)
{
    int disc_num = sim_InitPar0(theSim);

    if(disc_num == 0)
        cell_init_ntdisc_thompson(c, r, phi, z, theSim);
    else
        printf("ERROR: cell_init_ntdisc given bad option.\n");
}

void cell_single_init_ntdisc(struct Cell *theCell, struct Sim *theSim,int i,int j,int k)
{
    double rm = sim_FacePos(theSim,i-1,R_DIR);
    double rp = sim_FacePos(theSim,i,R_DIR);
    double r = 0.5*(rm+rp);
    double zm = sim_FacePos(theSim,k-1,Z_DIR);
    double zp = sim_FacePos(theSim,k,Z_DIR);
    double z = 0.5*(zm+zp);
    double t = theCell->tiph-.5*theCell->dphi;

    cell_init_ntdisc_calc(theCell, r, t, z, theSim);
}

void cell_init_ntdisc(struct Cell ***theCells,struct Sim *theSim,struct MPIsetup * theMPIsetup)
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
             
                cell_init_ntdisc_calc(&(theCells[k][i][j]), r, t, z, theSim);
                
                if(PRINTTOOMUCH)
                {
                    printf("(%d,%d,%d) = (%.12lg, %.12lg, %.12lg): (%.12lg, %.12lg, %.12lg, %.12lg, %.12lg)\n", i,j,k,r,t,z,theCells[k][i][j].prim[RHO],theCells[k][i][j].prim[URR],theCells[k][i][j].prim[UPP],theCells[k][i][j].prim[UZZ],theCells[k][i][j].prim[PPP]);
                }
            }
        }
    }
}
