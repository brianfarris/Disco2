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

  double omega = 0.0;
  double rm = sim_FacePos(theSim,i-1,R_DIR);
  double rp = sim_FacePos(theSim,i,R_DIR);
  double r_coord = 0.5*(rm+rp);
  double phi = theCell->tiph - .5*theCell->dphi;

  double x0 = 0.3;
  double x = r_coord*cos(phi)-x0;
  double y = r_coord*sin(phi);

  double Pp  = 1.0;
  double rho = 1.0;
  double B02 = 0.3*1e-3;

  double R = 0.25;

  double r = sqrt(x*x+y*y);

  double X = r/R;
  double B2 = 0.0; 
  double dP = -B02*((3./4.)*X*X*X*X - (5./3.)*X*X*X + X*X - (1./12.) );

  double Az = 0.0;

  if( X < 1.0 ){
    Pp += dP;
    B2 = B02*X*X*pow(1.-X,2.);
    Az = R*sqrt(B02)*( X*X*( .5 - X/3. ) -1./6. );
  }
  double Bx = -sqrt(B2)*y/r;
  double By =  sqrt(B2)*x/r; 
  theCell->prim[RHO] = rho;
  theCell->prim[PPP] = Pp;
  theCell->prim[URR] = 0.0;
  theCell->prim[UPP] = omega;
  theCell->prim[UZZ] = 0.0;
  theCell->prim[BRR] =  Bx*cos(phi) + By*sin(phi);
  theCell->prim[BPP] = -Bx*sin(phi) + By*cos(phi);
  theCell->prim[BZZ] = 0.0;
  theCell->prim[PSI] = 0.0;
  theCell->prim[ARR] = 0.0;
  theCell->prim[APP] = 0.0;
  theCell->prim[AZZ] = Az;
  theCell->prim[PHI] = 0.0;

  theCell->divB = 0.0;
  theCell->GradPsi[0] = 0.0;
  theCell->GradPsi[1] = 0.0;
  theCell->GradPsi[2] = 0.0;
}

void cell_init_fieldloop(struct Cell ***theCells,struct Sim *theSim,struct MPIsetup * theMPIsetup) {

  double R = 0.25;
  double B02 = 0.3*1e-3;
  double rho = 1.0;
  double P0 = 1.0;
  double x0 = 0.3;
  double omega = 0.0;

  int i, j,k;
  for (k = 0; k < sim_N(theSim,Z_DIR); k++) {
    for (i = 0; i < sim_N(theSim,R_DIR); i++) {
      double rm = sim_FacePos(theSim,i-1,R_DIR);
      double rp = sim_FacePos(theSim,i,R_DIR);
      double r_coord = 0.5*(rm+rp);
      double Pp = P0;// + .5*Rho0*Om*Om*r*r;

      for (j = 0; j < sim_N_p(theSim,i); j++) {
        double phi = theCells[k][i][j].tiph - .5*theCells[k][i][j].dphi;

        double x = r_coord*cos(phi)-x0;
        double y = r_coord*sin(phi);

        double r = sqrt(x*x+y*y);
        double X = r/R;
        double B2 = 0.0; 
        double dP = -B02*((3./4.)*X*X*X*X - (5./3.)*X*X*X + X*X - (1./12.) );

        double Az = 0.0;

        if( X < 1.0 ){
          Pp += dP;
          B2 = B02*X*X*pow(1.-X,2.);
          Az = R*sqrt(B02)*( X*X*( .5 - X/3. ) -1./6. );
        }
        double Bx = -sqrt(B2)*y/r;
        double By =  sqrt(B2)*x/r; 
        theCells[k][i][j].prim[RHO] = rho;
        theCells[k][i][j].prim[PPP] = Pp;
        theCells[k][i][j].prim[URR] = 0.0;
        theCells[k][i][j].prim[UPP] = omega;
        theCells[k][i][j].prim[UZZ] = 0.0;
        theCells[k][i][j].prim[BRR] =  Bx*cos(phi) + By*sin(phi);
        theCells[k][i][j].prim[BPP] = -Bx*sin(phi) + By*cos(phi);
        theCells[k][i][j].prim[BZZ] = 0.0;
        theCells[k][i][j].prim[PSI] = 0.0;
        theCells[k][i][j].prim[ARR] = 0.0;
        theCells[k][i][j].prim[APP] = 0.0;
        theCells[k][i][j].prim[AZZ] =  Az;
        theCells[k][i][j].prim[PHI] = 0.0; 

        theCells[k][i][j].divB = 0.0;
        theCells[k][i][j].GradPsi[0] = 0.0;
        theCells[k][i][j].GradPsi[1] = 0.0;
        theCells[k][i][j].GradPsi[2] = 0.0;
      }
    }
  }
}
