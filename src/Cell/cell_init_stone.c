#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Sim.h"
#include "../Headers/Face.h"
#include "../Headers/GravMass.h"
#include "../Headers/MPIsetup.h"
#include "../Headers/header.h"

void cell_single_init_stone(struct Cell *theCell, struct Sim *theSim,int i,int j,int k){
  double DISK_MACH = 20.;
  double GAMMALAW = sim_GAMMALAW(theSim);
  double rm =  sim_FacePos(theSim,i-1,R_DIR);
  double rp =  sim_FacePos(theSim,i,R_DIR);
  double r = .5*(rm+rp);
  double omega = 1./pow(r,1.5);
  double cs = 1.0/DISK_MACH*omega*r;
  double rho = 100.0;
  double Pp = cs*cs*rho/GAMMALAW;

  theCell->prim[RHO] = rho;
  theCell->prim[PPP] = Pp;
  theCell->prim[URR] = 0.0;
  theCell->prim[UPP] = 0.0;//omega;
  theCell->prim[UZZ] = 0.0;
  theCell->prim[BRR] = 0.0;
  theCell->prim[BPP] = 0.0;
  theCell->prim[BZZ] = 0.0;
  theCell->prim[PSI] = 0.0;
  theCell->wiph = omega*r;
  theCell->divB = 0.0;
  theCell->GradPsi[0] = 0.0;
  theCell->GradPsi[1] = 0.0;
  theCell->GradPsi[2] = 0.0;
  printf("you shouldn't need to call this\n");
  exit(1);
}

void cell_init_stone(struct Cell ***theCells,struct Sim *theSim,struct MPIsetup * theMPIsetup) {

  double rho0 = 100.0;
  double DISK_MACH = 20.;
  double GAMMALAW = sim_GAMMALAW(theSim);
  double H0     = 0.2;
  //double A_N    = H0*sqrt(rho0)*sqrt(15./16.)/(2.*M_PI);
  double A_N    = 0.01;

  srand(666 + mpisetup_MyProc(theMPIsetup));
  int i, j,k;
  for (k = 0; k < sim_N(theSim,Z_DIR); k++) {
    for (i = 0; i < sim_N(theSim,R_DIR); i++) {
      double rm = sim_FacePos(theSim,i-1,R_DIR);
      double rp = sim_FacePos(theSim,i,R_DIR);
      double r = 0.5*(rm+rp);
      double omega = 1./pow(r,1.5);
      double cs = 1.0/DISK_MACH;
      double rho = rho0;
      double Pp = cs*cs*rho/GAMMALAW;
      double IR = 1.0;
      if( r<1.2 || r>3.8 ) IR = 0.0;
      if( r<1.4 ) IR *= (r-1.2)/.2;
      if( r>3.6 ) IR *= (3.8-r)/.2;
      double Bz = A_N*IR*omega;
      if (BzZ==1) Bz *= sin(2*M_PI*(r-2.)/.2);

      double Ap = -A_N*IR*pow(r,-1.5)*cos(2*M_PI * (r-2.)/H0);

      //printf("r: %e, Ap: %e\n",r,Ap);

      for (j = 0; j < sim_N_p(theSim,i); j++) {
        double delta = .02*( (double)rand()/(double)RAND_MAX - .5 );

        theCells[k][i][j].prim[RHO] = rho;
        theCells[k][i][j].prim[PPP] = Pp;
        theCells[k][i][j].prim[URR] = 0.0;
        theCells[k][i][j].prim[UPP] = omega*delta;
        theCells[k][i][j].prim[UZZ] = delta;
        theCells[k][i][j].prim[BRR] = 0.0;
        theCells[k][i][j].prim[BPP] = 0.0;
        theCells[k][i][j].prim[BZZ] = 0.0;//Bz;
        theCells[k][i][j].prim[PSI] = 0.0;
        theCells[k][i][j].prim[ARR] = 0.0;
        theCells[k][i][j].prim[APP] = Ap;
        theCells[k][i][j].prim[AZZ] = 0.0;
        theCells[k][i][j].wiph = omega*r;
        theCells[k][i][j].divB = 0.0;
        theCells[k][i][j].GradPsi[0] = 0.0;
        theCells[k][i][j].GradPsi[1] = 0.0;
        theCells[k][i][j].GradPsi[2] = 0.0;
      }
    }
  }

}
