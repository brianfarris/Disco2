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

  double rho = pow(r,-3./5.);
  double   P = 1.e-2*pow(r,-3./2.);

  double o2  = 1./r/r/r;// - 3./2.*P/rho/r/r;
  double omega = sqrt(o2);
  double vr     = 0.0;//-1.5*alpha*Gam*(P/rho)*sqrt(r);


  theCell->prim[RHO] = rho;//*exp(-pow(r/20.,-2.));
  theCell->prim[PPP] = P;//*exp(-pow(r/20.,-2.));
  theCell->prim[URR] = vr;
  theCell->prim[UPP] = omega-1./pow(r,1.5);
  theCell->prim[UZZ] = 0.0;
  theCell->wiph = pow(r,-.5);
  theCell->drOm = -1.5*pow(r,-2.5);
  theCell->divB = 0.0;
  theCell->GradPsi[0] = 0.0;
  theCell->GradPsi[1] = 0.0;
  theCell->GradPsi[2] = 0.0;

}

void cell_init_middle(struct Cell ***theCells,struct Sim *theSim,struct MPIsetup * theMPIsetup) {

  double alpha = sim_EXPLICIT_VISCOSITY(theSim);
  double Gam = sim_GAMMALAW(theSim);

  int i, j,k;
  for (k = 0; k < sim_N(theSim,Z_DIR); k++) {
    for (i = 0; i < sim_N(theSim,R_DIR); i++) {
      double rm = sim_FacePos(theSim,i-1,R_DIR);
      double rp = sim_FacePos(theSim,i,R_DIR);
      double r = 0.5*(rm+rp);

      for (j = 0; j < sim_N_p(theSim,i); j++) {


        double rho = pow(r,-3./5.);//*exp(-pow(r/20.,-2.));
        double   P = 1.e-2*pow(r,-3./2.);//*exp(-pow(r/20.,-2.));
        double o2  = 1./r/r/r;// - 3./2.*P/rho/r/r;
        double omega = sqrt(o2);
        double vr     = 0.0;//-1.5*alpha*Gam*(P/rho)*sqrt(r);

        theCells[k][i][j].prim[RHO] = rho;
        theCells[k][i][j].prim[PPP] = P;
        theCells[k][i][j].prim[URR] = vr;
        theCells[k][i][j].prim[UPP] = omega-sqrt(1.)/pow(r,1.5);
        theCells[k][i][j].prim[UZZ] = 0.0;
        theCells[k][i][j].wiph = pow(r,-.5);
        theCells[k][i][j].drOm = -1.5*pow(r,-2.5);
        theCells[k][i][j].divB = 0.0;
        theCells[k][i][j].GradPsi[0] = 0.0;
        theCells[k][i][j].GradPsi[1] = 0.0;
        theCells[k][i][j].GradPsi[2] = 0.0;
      }
    }
  }
}
