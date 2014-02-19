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
  double Pp  = 0.01;
  double v0  = 1.0;
  double t0  = 0.25+time_global;

  double rm = sim_FacePos(theSim,i-1,R_DIR);
  double rp = sim_FacePos(theSim,i,R_DIR);
  double r = 0.5*(rm+rp);
  double t = theCell->tiph-.5*theCell->dphi;
  double x  = r*cos(t)-3.;
  double y = r*sin(t);
  
  double nu;
  if (sim_EXPLICIT_VISCOSITY(theSim)>0.0){
    nu = sim_EXPLICIT_VISCOSITY(theSim);
  } else{
    nu = 0.05;
  }
  /*
  double nu;
  if (sim_VISC_CONST(theSim)==1){
    nu = sim_EXPLICIT_VISCOSITY(theSim);
  } else{
  double HoR = 0.1;
  //nu = sim_EXPLICIT_VISCOSITY(theSim)*HoR*HoR*pow(fabs((r*cos(t))),1.5);
  nu = sim_EXPLICIT_VISCOSITY(theSim)*sim_GAMMALAW(theSim)*Pp/rho*pow(fabs(r/3.*cos(t)),1.5);
  }
  */

  double vy = 0.0;
  if (fabs(y)<50.){
    vy = v0/sqrt(2.*M_PI*nu*t0)*exp(-x*x/(4.*nu*t0));// /sqrt(2.*M_PI*sim_EXPLICIT_VISCOSITY(theSim)*t0);
  }

  double vr    = vy*sin(t);
  double omega = vy*cos(t)/r;

  theCell->prim[RHO] = rho;
  theCell->prim[PPP] = Pp;
  theCell->prim[URR] = vr;
  theCell->prim[UPP] = omega-sim_W_A(theSim,r)/r;
  theCell->prim[UZZ] = 0.0;
  theCell->wiph = 0.0;
  theCell->divB = 0.0;
  theCell->GradPsi[0] = 0.0;
  theCell->GradPsi[1] = 0.0;
  theCell->GradPsi[2] = 0.0;
}

void cell_init_shear(struct Cell ***theCells,struct Sim *theSim,struct MPIsetup * theMPIsetup) {

  double rho = 1.0;
  double Pp  = 0.01;
  double v0  = 1.0;
  double t0  = 0.25;

  int i, j,k;
  printf("sim_EXPLICIT_VISCOSITY(theSim): %e\n",sim_EXPLICIT_VISCOSITY(theSim));
  printf("sim_GAMMALAW(theSim): %e\n",sim_GAMMALAW(theSim));
  printf("sim_VISC_CONST(theSim): %d\n",sim_VISC_CONST(theSim));
  for (k = 0; k < sim_N(theSim,Z_DIR); k++) {
    for (i = 0; i < sim_N(theSim,R_DIR); i++) {
      double rm = sim_FacePos(theSim,i-1,R_DIR);
      double rp = sim_FacePos(theSim,i,R_DIR);
      double r = 0.5*(rm+rp);

      for (j = 0; j < sim_N_p(theSim,i); j++) {
        double t = theCells[k][i][j].tiph-.5*theCells[k][i][j].dphi;
        double x  = r*cos(t)-3.;
        double y = r*sin(t);
        /*
           double nu;
           if (sim_VISC_CONST(theSim)==1){
           nu = sim_EXPLICIT_VISCOSITY(theSim);
           } else{
           nu = sim_EXPLICIT_VISCOSITY(theSim)*sim_GAMMALAW(theSim)*Pp/rho*pow(fabs((r*cos(t))),1.5);
           }
           */
        double nu;
        if (sim_EXPLICIT_VISCOSITY(theSim)>0.0){
          nu = sim_EXPLICIT_VISCOSITY(theSim);
        } else{
          nu = 0.05;
        } 

        double vy=0.0;
        if (fabs(y)<50.){
          vy = v0/sqrt(2.*M_PI*nu*t0)*exp(-x*x/(4.*nu*t0)); // /sqrt(2.*M_PI*sim_EXPLICIT_VISCOSITY(theSim)*t0);
        } //else if (fabs(y)>=5. && fabs(y)<6.){
        // vy = v0/sqrt(2.*M_PI*nu*t0)*exp(-x*x/(4.*nu*t0)) * (6.-fabs(y));
        //}
        //printf("x: %e, x*x/(4.*sim_EXPLICIT_VISCOSITY(theSim)*t0): %e\n",x,x*x/(4.*sim_EXPLICIT_VISCOSITY(theSim)*t0));

        double vr    = vy*sin(t);
        double omega = vy*cos(t)/r;

        theCells[k][i][j].prim[RHO] = rho;
        theCells[k][i][j].prim[PPP] = Pp;
        theCells[k][i][j].prim[URR] = vr;
        theCells[k][i][j].prim[UPP] = omega-sim_W_A(theSim,r)/r;
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
