#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Sim.h"
#include "../Headers/Face.h"
#include "../Headers/GravMass.h"
#include "../Headers/header.h"


void cell_single_init_middle(struct Cell *theCell, struct Sim *theSim,int i,int j,int k){


  double rm = sim_FacePos(theSim,i-1,R_DIR);
  double rp = sim_FacePos(theSim,i,R_DIR);
  double r = 0.5*(rm+rp);

  double alpha = sim_EXPLICIT_VISCOSITY(theSim);
  double Gam = sim_GAMMALAW(theSim);
  
  double r0 = 2.5;
  int d  = 10;

  double rho = pow(r,-3./5.)*exp(-pow(r/r0,-d));
  double   P = PoRho_r1*pow(r,-3./2.)*exp(-pow(r/r0,-d));

  double r_drP_o_P = -1.5 + d*pow(r/r0,-d);

  if (rho<1.e-5) rho=1.e-5;
  if (P<1.e-8) P=1.e-8;
  double o2  = 1./r/r/r * pow( 1. + 3./16./r/r , 2) + r_drP_o_P * P/rho /r/r;
  double omega = sqrt(o2);
  double vr     = -3.*alpha/r/omega*P/rho*(2.+r_drP_o_P);


  theCell->prim[RHO] = rho;
  theCell->prim[PPP] = P;
  theCell->prim[URR] = vr;
  theCell->prim[UPP] = omega - sim_W_A(theSim,r)/r;
  theCell->prim[UZZ] = 0.0;
  theCell->wiph = 0.0;
  theCell->divB = 0.0;
  theCell->GradPsi[0] = 0.0;
  theCell->GradPsi[1] = 0.0;
  theCell->GradPsi[2] = 0.0;

}

void cell_init_middle(struct Cell ***theCells,struct Sim *theSim,struct MPIsetup * theMPIsetup) {

  double alpha = sim_EXPLICIT_VISCOSITY(theSim);
  double Gam = sim_GAMMALAW(theSim);

  double r0 = 2.5;
  int d = 10;
  int i, j,k;
  for (k = 0; k < sim_N(theSim,Z_DIR); k++) {
    for (i = 0; i < sim_N(theSim,R_DIR); i++) {
      double rm = sim_FacePos(theSim,i-1,R_DIR);
      double rp = sim_FacePos(theSim,i,R_DIR);
      double r = 0.5*(rm+rp);

      for (j = 0; j < sim_N_p(theSim,i); j++) {


        double rho = pow(r,-3./5.)*exp(-pow(r/r0,-d));
        double   P = PoRho_r1*pow(r,-3./2.)*exp(-pow(r/r0,-d));

        double r_drP_o_P = -1.5 + d*pow(r/r0,-d);

        if (rho<1.e-5) rho=1.e-5;
        if (P<1.e-7) P=1.e-7;

        double o2  = 1./r/r/r * pow( 1. + 3./16./r/r , 2) + r_drP_o_P * P/rho /r/r;
        double omega = sqrt(o2);
        double vr     = -3.*alpha/r/omega*P/rho*(2.+r_drP_o_P);

        theCells[k][i][j].prim[RHO] = rho;
        theCells[k][i][j].prim[PPP] = P;
        theCells[k][i][j].prim[URR] = vr;
        theCells[k][i][j].prim[UPP] = omega - sim_W_A(theSim,r)/r; 
        theCells[k][i][j].prim[UZZ] = 0.0;
        theCells[k][i][j].wiph = 0.0;
        theCells[k][i][j].divB = 0.0;
        theCells[k][i][j].GradPsi[0] = 0.0;
        theCells[k][i][j].GradPsi[1] = 0.0;
        theCells[k][i][j].GradPsi[2] = 0.0;
      }
    }
  }
}
