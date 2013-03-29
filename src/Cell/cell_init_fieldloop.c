#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Sim.h"
#include "../Headers/Face.h"
#include "../Headers/GravMass.h"
#include "../Headers/header.h"

void cell_single_init_fieldloop(struct Cell *theCell, struct Sim *theSim,int i,int j,int k){

  double Rl = 0.15;
  double B0 = 0.01;
  double Om = 10.0;
  double Rho0 = 1.0;
  double P0 = 1.0;
  double xinner = 1.0;
  double x0 = 0.3;
  double Vx = 0.0;

  double rm = sim_FacePos(theSim,i-1,R_DIR);
  double rp = sim_FacePos(theSim,i,R_DIR);
  double r = 0.5*(rm+rp);
  double omega = Om;
  double Pp = P0;// + .5*Rho0*Om*Om*r*r;

  double phi = theCell->tiph - .5*theCell->dphi;

  double x = r*cos(phi)-x0;
  double y = r*sin(phi);

  double rl = sqrt(x*x+y*y);
  double xx = M_PI*rl/Rl;

  double Bp = B0*pow(sin(xx-xinner),2.)*sqrt(2./M_PI*(xx-xinner));
  if( (xx-xinner > M_PI)||(xx<xinner) ) Bp = 0.0;

  double dP = B0*B0*( - (rl/Rl)*pow(sin(xx),4.) - (1./16./M_PI)*( 12.*xx - 8.*sin(2.*xx) + sin(4.*xx) ) );
  if ((xx-xinner > M_PI)||(xx<xinner) ) dP = 0.0;

  double Bx = -Bp*y/rl;
  double By =  Bp*x/rl;

  theCell->prim[RHO] = Rho0;
  theCell->prim[PPP] = Pp+dP;
  theCell->prim[URR] = Vx*cos(phi);
  theCell->prim[UPP] = omega-Vx*sin(phi);
  theCell->prim[UZZ] = 0.0;
  theCell->prim[BRR] =  Bx*cos(phi) + By*sin(phi);
  theCell->prim[BPP] = -Bx*sin(phi) + By*cos(phi);
  theCell->prim[BZZ] = 0.0;
  theCell->prim[PSI] = 0.0;

  theCell->divB = 0.0;
  theCell->GradPsi[0] = 0.0;
  theCell->GradPsi[1] = 0.0;
  theCell->GradPsi[2] = 0.0;
}

void cell_init_fieldloop(struct Cell ***theCells,struct Sim *theSim,struct MPIsetup * theMPIsetup) {

  double Rl = 0.15;
  double B0 = 0.01;
  double Om = 10.0;
  double Rho0 = 1.0;
  double P0 = 1.0;
  double xinner = 1.0;
  double x0 = 0.3;
  double Vx = 0.0;

  int i, j,k;
  for (k = 0; k < sim_N(theSim,Z_DIR); k++) {
    for (i = 0; i < sim_N(theSim,R_DIR); i++) {
      double rm = sim_FacePos(theSim,i-1,R_DIR);
      double rp = sim_FacePos(theSim,i,R_DIR);
      double r = 0.5*(rm+rp);
      double omega = Om;
      double Pp = P0;// + .5*Rho0*Om*Om*r*r;

      for (j = 0; j < sim_N_p(theSim,i); j++) {
        double phi = theCells[k][i][j].tiph - .5*theCells[k][i][j].dphi;

        double x = r*cos(phi)-x0;
        double y = r*sin(phi);

        double rl = sqrt(x*x+y*y);
        double xx = M_PI*rl/Rl;

        double Bp = B0*pow(sin(xx-xinner),2.)*sqrt(2./M_PI*(xx-xinner));
        if( (xx-xinner > M_PI)||(xx<xinner) ) Bp = 0.0;

        double dP = B0*B0*( - (rl/Rl)*pow(sin(xx),4.) - (1./16./M_PI)*( 12.*xx - 8.*sin(2.*xx) + sin(4.*xx) ) );
        if ((xx-xinner > M_PI)||(xx<xinner) ) dP = 0.0;

        double Bx = -Bp*y/rl;
        double By =  Bp*x/rl;

        theCells[k][i][j].prim[RHO] = Rho0;
        theCells[k][i][j].prim[PPP] = Pp+dP;
        theCells[k][i][j].prim[URR] = Vx*cos(phi);
        theCells[k][i][j].prim[UPP] = omega-Vx*sin(phi);
        theCells[k][i][j].prim[UZZ] = 0.0;
        theCells[k][i][j].prim[BRR] =  Bx*cos(phi) + By*sin(phi);
        theCells[k][i][j].prim[BPP] = -Bx*sin(phi) + By*cos(phi);
        theCells[k][i][j].prim[BZZ] = 0.0;
        theCells[k][i][j].prim[PSI] = 0.0;

        theCells[k][i][j].divB = 0.0;
        theCells[k][i][j].GradPsi[0] = 0.0;
        theCells[k][i][j].GradPsi[1] = 0.0;
        theCells[k][i][j].GradPsi[2] = 0.0;
      }
    }
  }
}
