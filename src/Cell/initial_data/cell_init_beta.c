#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../../Headers/Cell.h"
#include "../../Headers/Sim.h"
#include "../../Headers/Face.h"
#include "../../Headers/GravMass.h"
#include "../../Headers/header.h"


void cell_single_init_beta(struct Cell *theCell, struct Sim *theSim,struct GravMass * theGravMasses,int i,int j,int k){

    double r0 = 0.4;
    int d=2;
    double rm = sim_FacePos(theSim,i-1,R_DIR);
    double rp = sim_FacePos(theSim,i,R_DIR);
    double r = 0.5*(rm+rp);

    double Gamma = sim_GAMMALAW(theSim);
    double xi_g = 0.1;
    double xi_r_times_xi_g_pow_8 = 1.e-1;
    double xi_r =  xi_r_times_xi_g_pow_8*pow(xi_g,-8);

    double Sigma =  pow(r,-3./5.)*exp(-pow(r/r0,-d));
    if (Sigma<1.e-3) Sigma=1.e-3;
    double cs = 0.5*xi_r*pow(xi_g,8)*pow(r,-1.5)
        *(1+sqrt(1.+4.*xi_g*xi_g/xi_r/xi_r * pow(xi_g,-16) * pow(r,21./10.)))*exp(-.5*pow(r/r0,-d));
    double P = cs*cs*Sigma;
   // if (P<1.e-6) P=1.e-6;
    theCell->prim[RHO] = Sigma;
    theCell->prim[PPP] = P;
    theCell->prim[URR] = 0.0;
    theCell->prim[UPP] = 0.0;
    theCell->prim[UZZ] = 0.0;
    theCell->wiph = pow(r,-.5);
    //printf("ERROR. cell_single_init_beta isnt set up right now\n");
    //exit(0);
}

void cell_init_beta(struct Cell ***theCells,struct Sim *theSim,struct GravMass * theGravMasses,struct MPIsetup * theMPIsetup) {

    double r0 = 0.4;
    int d=2;
    double Gamma = sim_GAMMALAW(theSim);
    double xi_g = 0.1;
    double xi_r_times_xi_g_pow_8 = 1.e-1;
    double xi_r =  xi_r_times_xi_g_pow_8*pow(xi_g,-8);

    int i, j,k;
    for (k = 0; k < sim_N(theSim,Z_DIR); k++) {
        for (i = 0; i < sim_N(theSim,R_DIR); i++) {
            double rm = sim_FacePos(theSim,i-1,R_DIR);
            double rp = sim_FacePos(theSim,i,R_DIR);
            double r = 0.5*(rm+rp);

            for (j = 0; j < sim_N_p(theSim,i); j++) {

                double Sigma =  pow(r,-3./5.)*exp(-pow(r/r0,-d));
                if (Sigma<1.e-3) Sigma=1.e-3;
                double cs = 0.5*xi_r*pow(xi_g,8)*pow(r,-1.5)
                    *(1+sqrt(1.+4.*xi_g*xi_g/xi_r/xi_r * pow(xi_g,-16) * pow(r,21./10.)))*exp(-.5*pow(r/r0,-d));
                double P = cs*cs*Sigma;
                //if (P<1.e-6) P=1.e-6;


                //printf("%e %e\n",r,P);
                theCells[k][i][j].prim[RHO] = Sigma;
                theCells[k][i][j].prim[PPP] = P;
                theCells[k][i][j].prim[URR] = 0.0;
                theCells[k][i][j].prim[UPP] = 0.0;
                theCells[k][i][j].prim[UZZ] = 0.0;
                theCells[k][i][j].wiph = pow(r,-0.5);
               // printf("%e %e\n",r,P);
            }
        }
    }
   // exit(1);
}
