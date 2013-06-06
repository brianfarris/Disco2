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

void cell_single_init_torus(struct Cell *theCell, struct Sim *theSim,int i,int j,int k){
  printf("WARNING, YOU SHOULDNT BE HERE\n");
  exit(0);
}

void cell_init_torus(struct Cell ***theCells,struct Sim *theSim,struct MPIsetup * theMPIsetup) {
  double GAMMALAW = sim_GAMMALAW(theSim);
  double M = 1.0; 
  double r_inner = 0.7;
  double r_max = 1.0;
  double q = 2.;
  double s = 2.*q-q;
  double GammaFac = (GAMMALAW-1.)/GAMMALAW;
  double Phi_inner = -1.0/r_inner;
  double Phi_max = -1.0/r_max;
  double Rho_norm = pow(GammaFac*(Phi_inner-Phi_max+pow(r_inner,-s)/s-pow(1.0,-s)/s),1.0/(GAMMALAW-1.0));
  double PoRho_max = GammaFac*(Phi_inner-Phi_max+pow(r_inner,-s)/s-pow(1.0,-s)/s);

  int i, j,k;
  srand(666 + mpisetup_MyProc(theMPIsetup));
  for (k = 0; k < sim_N(theSim,Z_DIR); k++) {
    double z;
    if( sim_N_global(theSim,Z_DIR) != 1 ){
      double zp = sim_FacePos(theSim,k  ,Z_DIR);
      double zm = sim_FacePos(theSim,k-1,Z_DIR);
      z = 0.5*(zm+zp); 
    } else{
      z = 0.0;
    }
    for (i = 0; i < sim_N(theSim,R_DIR); i++) {
      double rm = sim_FacePos(theSim,i-1,R_DIR);
      double rp = sim_FacePos(theSim,i,R_DIR);
      double r = 0.5*(rm+rp);
      double R = sqrt( r*r + z*z );      
      double Phi = -M/R;
      double omega = pow(r,-q);
      if (r<r_inner){
        omega = pow(r_inner,-q)*r/r_inner;
      }
      double PoRho = GammaFac*(Phi_inner-Phi+pow(r_inner,-s)/s-pow(r,-s)/s);
      if (PoRho<1.0e-20){
        PoRho = 1.0e-20;
      }
      double rho = pow(PoRho,1.0/(GAMMALAW-1.0))/Rho_norm;
      double Pp = pow(PoRho,1.0/GammaFac)/Rho_norm;

      if (PoRho<1.0e-20){
        Pp = rho;
      }

      double magfac;
      if (rho>0.75){
        magfac = 0.05/PoRho_max;
      } else{
        magfac = 0.0/PoRho_max;
      }
      double Br = magfac*GammaFac*M/R/R*z/R;
      double Bp = 0.0;
      double Bz = magfac*(GammaFac*(-M/R/R*r/R+pow(r,-s-1.))+PoRho/r );


      for (j = 0; j < sim_N_p(theSim,i); j++) {
        double delta = .02*( (double)rand()/(double)RAND_MAX - .5 );
        theCells[k][i][j].prim[RHO] = rho;
        theCells[k][i][j].prim[PPP] = Pp;
        theCells[k][i][j].prim[URR] = 0.0;
        theCells[k][i][j].prim[UPP] = omega*(1.+delta);
        theCells[k][i][j].prim[UZZ] = 0.0;
        if (sim_runtype(theSim)==MHD){
        theCells[k][i][j].prim[BRR] = Br;
        theCells[k][i][j].prim[BPP] = Bp;
        theCells[k][i][j].prim[BZZ] = Bz;
        theCells[k][i][j].prim[PSI] = 0.0;
        theCells[k][i][j].divB = 0.0;
        theCells[k][i][j].GradPsi[0] = 0.0;
        theCells[k][i][j].GradPsi[1] = 0.0;
        theCells[k][i][j].GradPsi[2] = 0.0;
        }
        theCells[k][i][j].wiph = omega*r;
      }
    }
  }

}
