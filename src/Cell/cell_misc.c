#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Sim.h"
#include "../Headers/Face.h"
#include "../Headers/GravMass.h"
#include "../Headers/header.h"


void cell_clean_pi(struct Cell *** theCells,struct Sim *theSim){
  int NUM_Q = sim_NUM_Q(theSim);

  int i,j,k;
  for( k=0 ; k<sim_N(theSim,Z_DIR) ; ++k ){
    for( i=0 ; i<sim_N(theSim,R_DIR) ; ++i ){
      for( j=0 ; j<sim_N_p(theSim,i) ; ++j ){
        double phi = theCells[k][i][j].tiph;
        while( phi > 2.*M_PI ) phi -= 2.*M_PI;
        while( phi < 0.0 ) phi += 2.*M_PI;
        theCells[k][i][j].tiph = phi;
      }
    }
  }
}

void cell_copy(struct Cell ***theCells,struct Sim * theSim){
  int NUM_Q = sim_NUM_Q(theSim);
  int i,j,k,q;
  for( k=0 ; k<sim_N(theSim,Z_DIR) ; ++k ){
    for( i=0 ; i<sim_N(theSim,R_DIR) ; ++i ){
      for( j=0 ; j<sim_N_p(theSim,i) ; ++j ){
        for (q=0;q<NUM_Q;++q){
          theCells[k][i][j].RKcons[q] = theCells[k][i][j].cons[q];
        }
        theCells[k][i][j].RKtiph = theCells[k][i][j].tiph;
      }
    }
  }
}

void cell_adjust_RK_cons( struct Cell *** theCells , struct Sim * theSim, double RK ){
  int NUM_Q = sim_NUM_Q(theSim);

  int i,j,k,q;
  for( k=0 ; k<sim_N(theSim,Z_DIR) ; ++k ){
    for( i=0 ; i<sim_N(theSim,R_DIR) ; ++i ){
      for( j=0 ; j<sim_N_p(theSim,i) ; ++j ){
        struct Cell * c = &(theCells[k][i][j]);
        for( q=0 ; q<NUM_Q ; ++q ){
          c->cons[q] = (1.0-RK)*c->cons[q] + RK*c->RKcons[q];
        }
      }
    }
  }
}

void cell_update_phi( struct Cell *** theCells , struct Sim * theSim, double RK , double dt ){

  int i,j,k;
  for( k=0 ; k<sim_N(theSim,Z_DIR) ; ++k ){
    for( i=0 ; i<sim_N(theSim,R_DIR) ; ++i ){
      double rm = sim_FacePos(theSim,i-1,R_DIR);
      double rp = sim_FacePos(theSim,i,R_DIR);
      double r = .5*(rm+rp);
      for( j=0 ; j<sim_N_p(theSim,i) ; ++j ){
        struct Cell * c = &(theCells[k][i][j]);
        while( c->tiph - c->RKtiph >  M_PI ) c->RKtiph += 2.*M_PI;
        while( c->tiph - c->RKtiph < -M_PI ) c->RKtiph -= 2.*M_PI;
        c->tiph = (1.0-RK)*c->tiph + RK*c->RKtiph;
        double w = theCells[k][i][j].wiph;
        theCells[k][i][j].tiph += w*dt/r;
      }
    }
  }
}

void cell_update_dphi( struct Cell *** theCells,struct Sim * theSim ){

  int i,j,k;
  for( k=0 ; k<sim_N(theSim,Z_DIR) ; ++k ){
    for( i=0 ; i<sim_N(theSim,R_DIR) ; ++i ){
      for( j=0 ; j<sim_N_p(theSim,i) ; ++j ){
        int jm = j-1;
        if( jm==-1 ) jm = sim_N_p(theSim,i)-1;
        double dphi = theCells[k][i][j].tiph - theCells[k][i][jm].tiph;
        while( dphi > 2.*M_PI ) dphi -= 2.*M_PI;
        while( dphi < 0.0 ) dphi += 2.*M_PI;
        theCells[k][i][j].dphi = dphi;
      }
    }
  }
}

void cell_print(struct Cell *** theCells,int n){
  printf("%d theCells[2][5][4]: %e\n",n,theCells[2][5][4].wiph);
}
