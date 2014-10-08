#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Sim.h"
#include "../Headers/Face.h"
#include "../Headers/GravMass.h"
#include "../Headers/header.h"


void cell_set_w(struct Cell ***theCells,struct Sim *theSim){
  int i,j,k;
  if (sim_W_NUMERIC_TYPE(theSim) == C_WCELL){
    printf("THIS OPTION SHOULD NOT BE USED RIGHT NOW. WE NEED TO KNOW THE ANALYTIC RADIAL DERIV OF W\n");
    exit(1);
    for( k=0 ; k<sim_N(theSim,Z_DIR) ; ++k ){
      for( i=0 ; i<sim_N(theSim,R_DIR) ; ++i ){
        double rp = sim_FacePos(theSim,i,R_DIR);
        double rm = sim_FacePos(theSim,i-1,R_DIR);
        double r = 0.5*(rm+rp);
        for( j=0 ; j<sim_N_p(theSim,i) ; ++j ){
          int jp = j+1;
          double w  = theCells[k][i][j ].prim[UPP];
          double wp = theCells[k][i][jp].prim[UPP];
          theCells[k][i][j].wiph = .5*r*(w+wp);
        }
      }
    }

  } else if( sim_W_NUMERIC_TYPE(theSim) == C_RIGID ) {
    //printf("THIS OPTION SHOULD NOT BE USED RIGHT NOW. WE NEED TO KNOW THE ANALYTIC RADIAL DERIV OF W\n");
    //exit(1);
    for( k=0 ; k<sim_N(theSim,Z_DIR) ; ++k ){
      for( i=0 ; i<sim_N(theSim,R_DIR) ; ++i ){
        double w=0.0;
        double rp = sim_FacePos(theSim,i,R_DIR);
        double rm = sim_FacePos(theSim,i-1,R_DIR);
        double r = 0.5*(rm+rp);
        double Mring = 0.0;
        for( j=0 ; j<sim_N_p(theSim,i) ; ++j ){
          double vp = r*theCells[k][i][j].prim[UPP];
          double m = theCells[k][i][j].cons[DDD];
          Mring += m;
          w += vp*m;
        }
        for( j=0 ; j<sim_N_p(theSim,i) ; ++j ){
          theCells[k][i][j].wiph = w/Mring;   
        }
        //printf("%e %e\n",r,w/Mring);
      }
      //exit(1);
    }

  } else if (sim_W_NUMERIC_TYPE(theSim) == C_KEPLER ) {
    for( k=0 ; k<sim_N(theSim,Z_DIR) ; ++k ){
      for( i=0 ; i<sim_N(theSim,R_DIR) ; ++i ){
        double rp = sim_FacePos(theSim,i,R_DIR);
        double rm = sim_FacePos(theSim,i-1,R_DIR);
        double r = 0.5*(rm+rp);
        for( j=0 ; j<sim_N_p(theSim,i) ; ++j ){
          theCells[k][i][j].wiph = pow(r,-0.5);
        }
      }
    }
  } else if (sim_W_NUMERIC_TYPE(theSim) == C_OMEGA20) {
    for( k=0 ; k<sim_N(theSim,Z_DIR) ; ++k ){
      for( i=0 ; i<sim_N(theSim,R_DIR) ; ++i ){
        double rp = sim_FacePos(theSim,i,R_DIR);
        double rm = sim_FacePos(theSim,i-1,R_DIR);
        double r = 0.5*(rm+rp);
        for( j=0 ; j<sim_N_p(theSim,i) ; ++j ){
          theCells[k][i][j].wiph = r*20.;
        }
      }
    }
  } else if (sim_W_NUMERIC_TYPE(theSim) == C_MILOS) {
    for( k=0 ; k<sim_N(theSim,Z_DIR) ; ++k ){
      for( i=0 ; i<sim_N(theSim,R_DIR) ; ++i ){
        double rp = sim_FacePos(theSim,i,R_DIR);
        double rm = sim_FacePos(theSim,i-1,R_DIR);
        double r = 0.5*(rm+rp);
        for( j=0 ; j<sim_N_p(theSim,i) ; ++j ){
          theCells[k][i][j].wiph = (1.-exp(-pow(r,1.5)))/sqrt(r);
        }
      }
    }
  } else if (sim_W_NUMERIC_TYPE(theSim) == C_POWERLAW){
    for( k=0 ; k<sim_N(theSim,Z_DIR) ; ++k ){
      for( i=0 ; i<sim_N(theSim,R_DIR) ; ++i ){
        double rp = sim_FacePos(theSim,i,R_DIR);
        double rm = sim_FacePos(theSim,i-1,R_DIR);
        double r = 0.5*(rm+rp);
        for( j=0 ; j<sim_N_p(theSim,i) ; ++j ){
          theCells[k][i][j].wiph = pow(r,-1.0);
        }
      }
    }

  }

}
