#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../../Headers/Cell.h"
#include "../../Headers/Sim.h"
#include "../../Headers/Face.h"
#include "../../Headers/GravMass.h"
#include "../../Headers/MPIsetup.h"
#include "../../Headers/header.h"

void cell_single_init_flock(struct Cell *theCell, struct Sim *theSim,int i,int j,int k){
  double DISK_MACH = 10.;
  double GAMMALAW = sim_GAMMALAW(theSim);
  double rm =  sim_FacePos(theSim,i-1,R_DIR);
  double rp =  sim_FacePos(theSim,i,R_DIR);
  double r = .5*(rm+rp);
  double omega = 1./pow(r,1.5);
  double cs = 1.0/DISK_MACH;//*omega*r;
  double rho = 1.0;
  double Pp = cs*cs*rho/GAMMALAW;

  theCell->prim[RHO] = 1.0;
  theCell->prim[PPP] = Pp;
  theCell->prim[URR] = 0.0;
  theCell->prim[UPP] = 0.0;//omega;
  theCell->prim[UZZ] = 0.0;
  theCell->prim[BRR] = 0.0;
  theCell->prim[BPP] = 0.0;
  theCell->prim[BZZ] = 0.0;
  theCell->prim[PSI] = 0.0;
  theCell->wiph = omega*r;
}

void cell_init_flock(struct Cell ***theCells,struct Sim *theSim,struct MPIsetup * theMPIsetup) {

  double DISK_MACH = 10.;
  double GAMMALAW = sim_GAMMALAW(theSim);
  
  srand(666 + mpisetup_MyProc(theMPIsetup));
  double rho0 = 1.0;
  int i, j,k;
  for (k = 0; k < sim_N(theSim,Z_DIR); k++) {
    for (i = 0; i < sim_N(theSim,R_DIR); i++) {
      double rm = sim_FacePos(theSim,i-1,R_DIR);
      double rp = sim_FacePos(theSim,i,R_DIR);
      double r = 0.5*(rm+rp);
      
      double omega = 1./pow(r,1.5);
      double cs = 1.0/DISK_MACH;//*omega*r;
      double rho = rho0;
      double Pp = cs*cs*rho/GAMMALAW;

      double Bz;
      double delta=0.0;
      if ((r>2.)&&(r<3.)){
        double n=4.0;
        Bz = .054485/n;
        delta = .001*( (double)rand()/(double)RAND_MAX - .5 );
      } else{
        Bz = 0.0;
      }

      for (j = 0; j < sim_N_p(theSim,i); j++) {
        theCells[k][i][j].prim[RHO] = rho;
        theCells[k][i][j].prim[PPP] = Pp;
        theCells[k][i][j].prim[URR] = 0.0;
        theCells[k][i][j].prim[UPP] = omega*(/*1.+*/delta);
        theCells[k][i][j].prim[UZZ] = delta;
        theCells[k][i][j].prim[BRR] = 0.0;
        theCells[k][i][j].prim[BPP] = 0.0;
        theCells[k][i][j].prim[BZZ] = Bz;
        theCells[k][i][j].prim[PSI] = 0.0;
        theCells[k][i][j].wiph = omega*r;
      }
    }
  }

}
