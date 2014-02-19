#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Sim.h"
#include "../Headers/Face.h"
#include "../Headers/GravMass.h"
#include "../Headers/header.h"

// Bondi Accretion in Schwarzschild coordinates with a pulse.
void cell_single_init_gbondi2(struct Cell *theCell, struct Sim *theSim,int i,int j,int k)
{
    double rho, Pp, vr, u, T;
    double GAMMALAW = sim_GAMMALAW(theSim);
    double RG = sim_GravRadius(theSim);
    int dim = sim_InitPar0(theSim)+2;
    double rs = sim_InitPar1(theSim);
    double Mdot = sim_InitPar2(theSim);
    double EPS = 1.0e-10;

    double rm = sim_FacePos(theSim,i-1,R_DIR);
    double rp = sim_FacePos(theSim,i,R_DIR);
    double r = 0.5*(rm+rp);
    double zm = sim_FacePos(theSim,k-1,Z_DIR);
    double zp = sim_FacePos(theSim,k,Z_DIR);
    double z = 0.5*(zm+zp);
    double t = theCell->tiph-.5*theCell->dphi;

    double n = 1.0/(GAMMALAW-1.0);
    double us, Ts, K;
    double C1, C2;

    if(dim==2)
    {
        us = sqrt(0.5*RG/rs);
        Ts = n*us*us / ((1+n)*(1-(n+1)*us*us));
        C1 = us * pow(Ts,n) *rs;
        K = pow(2*M_PI*C1/Mdot, GAMMALAW-1);
    }
    else if(dim==3)
    {
        us = sqrt(0.25*RG/rs);
        Ts = n*us*us / ((1+n)*(1-(n+3)*us*us));
        C1 = us * pow(Ts,n) *rs*rs;
        K = pow(4*M_PI*C1/Mdot, GAMMALAW-1);
    }
    else
        printf("ERROR: Bad value in cell_init_gbondi.InitPar0 = %d\n", dim-2);

    C2 = (1.0+(1.0+n)*Ts)*(1+(1+n)*Ts)*(1-RG/rs+us*us);

    double T0 = (sqrt(C2)-1) / (1+n);

    double c1 = 1.0 + n;
    double c2 = 1.0 - RG/r;
    double c3;
    if(dim==2)
        c3 = C1*C1/(r*r);
    else if (dim==3)
        c3 = C1*C1/(r*r*r*r);

    double t0 = -10;
    double t1 = T0;

    int count = 0;
    while(fabs((t1-t0)/t0) > EPS)
    {
        t0 = t1;
        double f1 = (1.0+c1*t0)*(1.0+c1*t0);
        double df1 = 2.0*c1*(1.0+c1*t0);
        double f2 = c2 + c3*pow(t0,-2.0*n);
        double df2 = -2.0*n*c3*pow(t0,-2.0*n-1.0);
        double f = f1*f2 - C2;
        double df = df1*f2 + df2*f1;

        if(r < rs && df > 0.0)
            t1 = 0.5*t0;
        else if(r >= rs && df < 0.0)
            t1 = 2.0*t0;
        else
            t1 -= f/df;
        
        if(t1 < 0)
            t1 = 0.5*t0;

        count++;
        if(count > 10000)
        {
            printf("Bondi BC failed to converge. (r=%lg)\n",r);
            break;
        }
    }
    T = t1;
    if (dim == 2)
        u = C1/(pow(T,n)*r);
    else if (dim == 3)
        u = C1/(pow(T,n)*r*r);
    
    double rho0 = pow(T0/K,n);
    double Pp0 = rho0 * T0;

    rho = pow(T/K, n);
    Pp = rho * T;
    vr = - sqrt((1.0-RG/r) / (1+u*u/(1.0-RG/r))) * u;

    theCell->prim[RHO] = rho;
    theCell->prim[PPP] = Pp;
    theCell->prim[URR] = vr;
    theCell->prim[UPP] = 0.0;
    theCell->prim[UZZ] = 0.0;
    theCell->divB = 0.0;
    theCell->GradPsi[0] = 0.0;
    theCell->GradPsi[1] = 0.0;

  //TODO: Not sure what this is for.  Ask someone if important.
  //if(sim_NUM_C(theSim)<sim_NUM_Q(theSim)) theCell->prim[sim_NUM_C(theSim)] = Qq;
}

void cell_init_gbondi2(struct Cell ***theCells,struct Sim *theSim,struct MPIsetup * theMPIsetup)
{

    double rho, Pp, vr, u, T;
    double GAMMALAW = sim_GAMMALAW(theSim);
    double RG = sim_GravRadius(theSim);
    int dim = sim_InitPar0(theSim)+2;
    double rs = sim_InitPar1(theSim);
    double Mdot = sim_InitPar2(theSim);
    double EPS = 1.0e-10;

    double n = 1.0/(GAMMALAW-1.0);

    double us, Ts, K;
    double C1, C2;

    if(dim == 2)
    {
        us = sqrt(0.5*RG/rs);
        Ts = n*us*us / ((1+n)*(1-(n+1)*us*us));
        C1 = us * pow(Ts,n) *rs;
        K = pow(2.0*M_PI*C1/Mdot, GAMMALAW-1);
    }
    else if (dim == 3)
    {
        us = sqrt(0.25*RG/rs);
        Ts = n*us*us / ((1+n)*(1-(n+3)*us*us));
        C1 = us * pow(Ts,n) *rs*rs;
        K = pow(4.0*M_PI*C1/Mdot, GAMMALAW-1);
    }
    else
        printf("ERROR: Bad value in cell_init_gbondi. InitPar0 = %d\n", dim-2);

    C2 = (1.0+(1.0+n)*Ts)*(1+(1+n)*Ts)*(1-RG/rs+us*us);
            
    double T0 = (sqrt(C2)-1) / (1+n);
    T = T0;

    double rho0 = pow(T0/K,n);
    double Pp0 = rho0 * T0;
    printf("Setting up Bondi Problem: dim=%d, RG=%lg, rs=%lg, Mdot=%lg\n",dim,RG,rs,Mdot);
    printf("us: %lg, Ts: %lg, K: %lg, C1: %lg\n", us, Ts, K, C1);
    printf("rho0: %lg, P0: %lg\n", rho0, Pp0);

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

            double c1 = 1.0 + n;
            double c2 = 1.0 - RG/r;
            double c3;
            if(dim == 2)
                c3 = C1*C1/(r*r);
            else if (dim == 3)
                c3 = C1*C1/(r*r*r*r);

            double t0 = -10;
            double t1;
            if(r < rs)
                t1 = 2*Ts;
            else
                t1 = 0.5*Ts;

            double f1,f2,df1,df2,f,df;
            int count = 0;
            while(fabs((t1-t0)/t0) > EPS)
            {
                t0 = t1;

                f1 = (1.0+c1*t0)*(1.0+c1*t0);
                df1 = 2.0*c1*(1.0+c1*t0);
                f2 = c2 + c3*pow(t0,-2.0*n);
                df2 = -2.0*n*c3*pow(t0,-2.0*n-1.0);
                f = f1*f2 - C2;
                df = df1*f2 + df2*f1;
                
                if(r < rs && df > 0.0)
                    t1 = 0.5*t0;
                else if(r >= rs && df < 0.0)
                    t1 = 2.0*t0;
                else
                    t1 -= f/df;
                
                if(t1 < 0)
                    t1 = 0.5*t0;

                count++;
                if(count > 100)
                {
                    printf("Bondi BC failed to converge. (r=%lg, f=%lg, df=%lg)\n",r,f,df);
                    break;
                }
            }
            T = t1;
            if(dim == 2)
                u = C1/(pow(T,n)*r);
            else if (dim == 3)
                u = C1/(pow(T,n)*r*r);

            for (j = 0; j < sim_N_p(theSim,i); j++) 
            {
                double t = theCells[k][i][j].tiph-.5*theCells[k][i][j].dphi;
             
            	rho = pow(T/K, n) + 3*rho0*exp(-((r-2.5*rs)*(r-2.5*rs)+r*r*(t-0.5*M_PI)*(t-0.5*M_PI))/(0.01*rs*rs));
            	Pp = rho * T;
            	vr = - sqrt((1.0-RG/r) / (1+u*u/(1.0-RG/r))) * u;

                theCells[k][i][j].prim[RHO] = rho;
                theCells[k][i][j].prim[PPP] = Pp;
                theCells[k][i][j].prim[URR] = vr;
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
