#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Sim.h"
#include "../Headers/header.h"


void cell_setT( struct Cell *** theCells ,struct Sim * theSim ){
  int i,j,k;
  double HoR = 0.05;
  double T0 = 0.005;
  double r_a = sqrt(1./log(T0*sim_GAMMALAW(theSim)/(HoR*HoR)));
  for( k=0 ; k<sim_N(theSim,Z_DIR) ; ++k ){
    double zm = sim_FacePos(theSim,k-1,Z_DIR);
    double zp = sim_FacePos(theSim,k,Z_DIR);
    for( i=0 ; i<sim_N(theSim,R_DIR) ; ++i ){
      double rm = sim_FacePos(theSim,i-1,R_DIR);
      double rp = sim_FacePos(theSim,i,R_DIR);
      double r = 0.5*(rm+rp);
      for( j=0 ; j<sim_N_p(theSim,i) ; ++j ){
        struct Cell * c = &(theCells[k][i][j]);
        double HoR = 0.05;
	double Omega_K;
	double cs;
        double PoRho;
        if (r>1.){
		Omega_K = pow(r,-1.5);
        	cs = HoR*Omega_K*r;
                PoRho = cs*cs/sim_GAMMALAW(theSim);
        } else{
		Omega_K = 1.;
		cs = HoR;
                PoRho = T0 * exp(-r*r/(r_a*r_a)); 
	}
        double rho = c->prim[RHO];
        c->prim[PPP] = PoRho * rho;
      }
    }
  }
}
