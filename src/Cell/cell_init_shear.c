#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Sim.h"
#include "../Headers/Face.h"
#include "../Headers/GravMass.h"
#include "../Headers/header.h"

void cell_single_init_shear(struct Cell ***theCells, struct Sim *theSim,int i,int j,int k){
  double rho = 1.0;
  double Pp  = 1.0;
  double v0  = 0.1;
  double t0  = 0.2;

  double nu = sim_EXPLICIT_VISCOSITY(theSim);
  double rm = sim_FacePos(theSim,i-1,R_DIR);
  double rp = sim_FacePos(theSim,i,R_DIR);
  double r = 0.5*(rm+rp);
  double t = theCells[k][i][j].tiph-.5*theCells[k][i][j].dphi;
  double x  = r*cos(t)-1.;

  double vy = v0*exp(-x*x/(2.*nu*t0))/sqrt(2.*M_PI*nu*t0);

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

void cell_init_shear(struct Cell ***theCells,struct Sim *theSim,struct MPIsetup * theMPIsetup) {

  double rho = 1.0;
  double Pp  = 1.0;
  double v0  = 0.1;
  double t0  = 0.2;

  double nu = sim_EXPLICIT_VISCOSITY(theSim);
  int i, j,k;

  for (k = 0; k < sim_N_z(theSim); k++) {
    for (i = 0; i < sim_N_r(theSim); i++) {
      double rm = sim_FacePos(theSim,i-1,R_DIR);
      double rp = sim_FacePos(theSim,i,R_DIR);
      double r = 0.5*(rm+rp);
      for (j = 0; j < sim_N_p(theSim,i); j++) {
        double t = theCells[k][i][j].tiph-.5*theCells[k][i][j].dphi;
        double x  = r*cos(t)-1.;

        double vy = v0*exp(-x*x/(2.*nu*t0))/sqrt(2.*M_PI*nu*t0);

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
