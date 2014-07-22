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
    double Mdot = sim_InitPar1(theSim);
    double GAM = sim_GAMMALAW(theSim);
    double M = sim_GravM(theSim);
    double alpha = sim_AlphaVisc(theSim);

    double rho, vr, vp, Pp;

    double eps = (5.0/3.0 - GAM) / (GAM - 1.0); // eps' parameter from Narayan & Yi
    double g = sqrt(1.0+18.0*alpha*alpha*GAM/((5+2*eps)*(5+2*eps))) - 1.0;

    //TODO: Make relativistic?  Confirm formulae.
    vr = -(5+2*eps)*g/(3*alpha*sqrt(GAM)) * sqrt(M/r);
    vp = sqrt(2*eps*(5+2*eps)*g/(9.0*alpha*alpha*GAM)) * sqrt(M/(r*r*r));

    double u0 = 1.0/sqrt(1.0-2*M/r-vr*vr/(1.0-2*M/r)-r*r*vp*vp);

    rho = -Mdot/(2*M_PI*r*u0*vr); //Vertically-integrated Density
    Pp = rho*(2*(5+2*eps)*g/(9*alpha*alpha*GAM))*(M/r);  //Vertically-integrated Pressure
 
    eps = (3.0 - GAM) / (GAM - 1.0); // eps' parameter from Narayan & Yi
    g = sqrt(1.0+2*18.0*18.0*alpha*alpha*GAM/((9+2*eps)*(9+2*eps))) - 1.0;

    //TODO: Make relativistic?  Confirm formulae.
    vr = -(9+2*eps)*g/(18*alpha*sqrt(GAM)) * sqrt(M/r);
    vp = sqrt(eps*(9+2*eps)*g/(2*81*alpha*alpha*GAM)) * sqrt(M/(r*r*r));
    u0 = 1.0/sqrt(1.0-2*M/r-vr*vr/(1.0-2*M/r)-r*r*vp*vp);

    rho = -Mdot/(2*M_PI*r*u0*vr); //Vertically-integrated Density
    Pp = rho * ((9+2*eps)*g/(54*alpha*alpha*GAM)) * (M/r);  //Vertically-integrated Pressure
    
    if(sim_Metric(theSim) == SCHWARZSCHILD_KS)
    {
        if(r < 10)
            vr = -2*M/(r+2*M);
        else
            vr = -Mdot/(2*M_PI*rho*r);

        if(r>2*M)
            vp = exp(-1.0/(r/M-2.0)) * sqrt(M/(r*r*r));
        else
            vp = 0.0;
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
