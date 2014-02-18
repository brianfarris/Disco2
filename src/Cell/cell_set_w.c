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
  if (sim_MOVE_CELLS(theSim) == C_WCELL){
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

  } else if( sim_MOVE_CELLS(theSim) == C_RIGID ) {
    printf("THIS OPTION SHOULD NOT BE USED RIGHT NOW. WE NEED TO KNOW THE ANALYTIC RADIAL DERIV OF W\n");
    exit(1);
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
      }
    }

  } else if (sim_MOVE_CELLS(theSim) == C_KEPLER ) {
    for( k=0 ; k<sim_N(theSim,Z_DIR) ; ++k ){
      for( i=0 ; i<sim_N(theSim,R_DIR) ; ++i ){
        double rp = sim_FacePos(theSim,i,R_DIR);
        double rm = sim_FacePos(theSim,i-1,R_DIR);
        double r = 0.5*(rm+rp);
        for( j=0 ; j<sim_N_p(theSim,i) ; ++j ){
          theCells[k][i][j].wiph = pow(r,-0.5);
          theCells[k][i][j].drOm = -1.5*pow(r,-2.5);
          //theCells[k][i][j].wiph = r/pow(r*r+sim_G_EPS(theSim)*sim_G_EPS(theSim),3./4.);
          //theCells[k][i][j].drOm = -1.5*r/pow(r*r+sim_G_EPS(theSim)*sim_G_EPS(theSim),7./4.);
        }
      }
    }
  } else if (sim_MOVE_CELLS(theSim) == C_OMEGA20) {
    for( k=0 ; k<sim_N(theSim,Z_DIR) ; ++k ){
      for( i=0 ; i<sim_N(theSim,R_DIR) ; ++i ){
        double rp = sim_FacePos(theSim,i,R_DIR);
        double rm = sim_FacePos(theSim,i-1,R_DIR);
        double r = 0.5*(rm+rp);
        for( j=0 ; j<sim_N_p(theSim,i) ; ++j ){
          theCells[k][i][j].wiph = r*20.;
          theCells[k][i][j].drOm = 20.;
        }
      }
    }
  } else if (sim_MOVE_CELLS(theSim) == C_MILOS) {
    for( k=0 ; k<sim_N(theSim,Z_DIR) ; ++k ){
      for( i=0 ; i<sim_N(theSim,R_DIR) ; ++i ){
        double rp = sim_FacePos(theSim,i,R_DIR);
        double rm = sim_FacePos(theSim,i-1,R_DIR);
        double r = 0.5*(rm+rp);
        for( j=0 ; j<sim_N_p(theSim,i) ; ++j ){
          /*
          if (r<0.1){
            theCells[k][i][j].wiph = pow(0.1,-0.5)*r;
            theCells[k][i][j].drOm = 0.0;
          } else{
            theCells[k][i][j].wiph = pow(r,-0.5);
            theCells[k][i][j].drOm = -1.5*pow(r,-2.5);
          }
          */
          theCells[k][i][j].wiph = (1.-exp(-pow(r,1.5)))/sqrt(r);
          theCells[k][i][j].drOm = 1.5*pow(r,-2.5)*(-1.+(1.+pow(r,1.5))*exp(-pow(r,1.5))); 
        }
      }
    }
  } else if (sim_MOVE_CELLS(theSim) == C_POWERLAW){
    for( k=0 ; k<sim_N(theSim,Z_DIR) ; ++k ){
      for( i=0 ; i<sim_N(theSim,R_DIR) ; ++i ){
        double rp = sim_FacePos(theSim,i,R_DIR);
        double rm = sim_FacePos(theSim,i-1,R_DIR);
        double r = 0.5*(rm+rp);
        for( j=0 ; j<sim_N_p(theSim,i) ; ++j ){
          theCells[k][i][j].wiph = pow(r,-1.0);
          theCells[k][i][j].drOm = -2.*pow(r,-3.0);
        }
      }
    }

  }

}
