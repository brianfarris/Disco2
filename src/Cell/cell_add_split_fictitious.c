#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Sim.h"
#include "../Headers/GravMass.h"
#include "../Headers/header.h"


void cell_add_split_fictitious( struct Cell *** theCells ,struct Sim * theSim,double dt){
  int i,j,k;
  for( k=0 ; k<sim_N(theSim,Z_DIR) ; ++k ){
    double zm = sim_FacePos(theSim,k-1,Z_DIR);
    double zp = sim_FacePos(theSim,k,Z_DIR);
    double dz = zp-zm;
    for( i=0 ; i<sim_N(theSim,R_DIR) ; ++i ){
      double rm = sim_FacePos(theSim,i-1,R_DIR);
      double rp = sim_FacePos(theSim,i,R_DIR);
      double r = 0.5*(rm+rp);
      for( j=0 ; j<sim_N_p(theSim,i) ; ++j ){
        struct Cell * c = &(theCells[k][i][j]);
        double dphi = c->dphi;
        double dV = dphi*.5*(rp*rp-rm*rm)*dz;
        double Omega = sim_W_A(theSim,r)/r;
        double dr_r2Om = 2.*sim_W_A(theSim,r) + r*r*sim_OM_A_DERIV(theSim,r);

        double kappa2 = 2.*Omega/r*dr_r2Om;
        if (kappa2>0.){
          double kappa = sqrt(kappa2);
          double l0 = c->cons[LLL];
          double rhov0 = c->cons[SRR];

          double l = l0*cos(kappa*dt) - c->cons[SRR]*dr_r2Om/kappa*sin(kappa*dt);
          double rhov = rhov0*cos(kappa*dt) + 2.0*Omega*c->cons[LLL]/r/kappa*sin(kappa*dt);

          c->cons[LLL] = l;
          c->cons[SRR] = rhov;
        }
      }
    }
  }
}
