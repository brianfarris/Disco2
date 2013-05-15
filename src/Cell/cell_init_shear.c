#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Sim.h"
#include "../Headers/Face.h"
#include "../Headers/GravMass.h"
#include "../Headers/header.h"

void cell_single_init_shear(struct Cell *theCell, struct Sim *theSim,int i,int j,int k){
  double rho = 1.0;
  double Pp  = 1.0;
  double v0  = 1.0;
  double t0  = 0.1+time_global;

  double rm = sim_FacePos(theSim,i-1,R_DIR);
  double rp = sim_FacePos(theSim,i,R_DIR);
  double r = 0.5*(rm+rp);
  double t = theCell->tiph-.5*theCell->dphi;
  double x  = r*cos(t)-4.;
  double nu;
  if (VISC_CONST==1){
    nu = sim_EXPLICIT_VISCOSITY(theSim);
  } else{
    nu = sim_EXPLICIT_VISCOSITY(theSim)*sim_GAMMALAW(theSim)*Pp/rho*pow(fabs(r/4.*cos(t)),1.5);
  }

  double vy = v0*exp(-x*x/(4.*sim_EXPLICIT_VISCOSITY(theSim)*t0))/sqrt(2.*M_PI*sim_EXPLICIT_VISCOSITY(theSim)*t0);

  double vr    = vy*sin(t);
  double omega = vy*cos(t)/r;

  theCell->prim[RHO] = rho;
  theCell->prim[PPP] = Pp;
  theCell->prim[URR] = vr;
  theCell->prim[UPP] = omega;
  theCell->prim[UZZ] = 0.0;
  theCell->wiph = 0.0;
  theCell->divB = 0.0;
  theCell->GradPsi[0] = 0.0;
  theCell->GradPsi[1] = 0.0;
  theCell->GradPsi[2] = 0.0;
}

void cell_init_shear(struct Cell ***theCells,struct Sim *theSim,struct MPIsetup * theMPIsetup) {

  double rho = 1.0;
  double Pp  = 1.0;
  double v0  = 1.0;
  double t0  = 0.1;

  int i, j,k;

  for (k = 0; k < sim_N(theSim,Z_DIR); k++) {
    for (i = 0; i < sim_N(theSim,R_DIR); i++) {
      double rm = sim_FacePos(theSim,i-1,R_DIR);
      double rp = sim_FacePos(theSim,i,R_DIR);
      double r = 0.5*(rm+rp);

     for (j = 0; j < sim_N_p(theSim,i); j++) {
       double t = theCells[k][i][j].tiph-.5*theCells[k][i][j].dphi;
       double x  = r*cos(t)-4.;

       double nu;
       if (VISC_CONST==1){
         nu = sim_EXPLICIT_VISCOSITY(theSim);
       } else{
         nu = sim_EXPLICIT_VISCOSITY(theSim)*sim_GAMMALAW(theSim)*Pp/rho*pow(fabs(r/4.*cos(t)),1.5);
       }

       double vy = v0*exp(-x*x/(4.*sim_EXPLICIT_VISCOSITY(theSim)*t0))/sqrt(2.*M_PI*sim_EXPLICIT_VISCOSITY(theSim)*t0);

       double vr    = vy*sin(t);
       double omega = vy*cos(t)/r;

       theCells[k][i][j].prim[RHO] = rho;
       theCells[k][i][j].prim[PPP] = Pp;
       theCells[k][i][j].prim[URR] = vr;
       theCells[k][i][j].prim[UPP] = omega;
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
