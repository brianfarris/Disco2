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

void cell_bc_damp( struct Cell *** theCells , struct Sim * theSim, double dt,void (*single_init_ptr)(struct Cell *,struct Sim *,int,int,int) ){
  int NUM_Q = sim_NUM_Q(theSim);

  double RMIN = sim_MIN(theSim,R_DIR);
  double RMAX = sim_MAX(theSim,R_DIR);
  double ZMIN = sim_MIN(theSim,Z_DIR);
  double ZMAX = sim_MAX(theSim,Z_DIR);

  struct Cell * initialCell = cell_single_create(theSim);
  double R0 = 1./sim_DAMP_TIME(theSim)/2./M_PI;

  int i,j,k,q;
  for( i=0 ; i<sim_N(theSim,R_DIR) ; ++i ){
    double rp = sim_FacePos(theSim,i,R_DIR);
    double rm = sim_FacePos(theSim,i-1,R_DIR);
    double r = .5*(rp+rm);

    if( r < sim_RDAMP_INNER(theSim) || r > sim_RDAMP_OUTER(theSim) ){
      double rate = R0*(r-sim_RDAMP_INNER(theSim))/(RMIN-sim_RDAMP_INNER(theSim));
      if( r > sim_RDAMP_OUTER(theSim) ) rate = R0*(r-sim_RDAMP_OUTER(theSim))/(RMAX-sim_RDAMP_OUTER(theSim));

      (*single_init_ptr)(initialCell,theSim,i,j,k);
      for( k=0 ; k<sim_N(theSim,Z_DIR) ; ++k ){
        for( j=0 ; j<sim_N_p(theSim,i) ; ++j ){
          for( q=0 ; q<NUM_Q ; ++q ){
            double dprim = initialCell->prim[q] - theCells[k][i][j].prim[q];
            theCells[k][i][j].prim[q] += (1.-exp(-rate*dt))*dprim;
          }
        }
      }
    }
  }
}

void cell_print(struct Cell *** theCells,int i,int j,int k){
  printf("at i: %d, j: %d, k: %d, DDD: %e\n",i,j,k,theCells[k][i][j].cons[DDD]);
  printf("at i: %d, j: %d, k: %d, TAU: %e\n",i,j,k,theCells[k][i][j].cons[TAU]);
  printf("at i: %d, j: %d, k: %d, SRR: %e\n",i,j,k,theCells[k][i][j].cons[SRR]);
  printf("at i: %d, j: %d, k: %d, LLL: %e\n",i,j,k,theCells[k][i][j].cons[LLL]);
  printf("at i: %d, j: %d, k: %d, SZZ: %e\n",i,j,k,theCells[k][i][j].cons[SZZ]);
}  
