#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Grid.h"
#include "../Headers/Face.h"
#include "../Headers/GravMass.h"
#include "../Headers/header.h"


void cell_clean_pi(struct Cell *** theCells,struct Grid *theGrid){
  int N_r_withghost = grid_N_r(theGrid)+grid_Nghost_rmin(theGrid)+grid_Nghost_rmax(theGrid);
  int N_z_withghost = grid_N_z(theGrid)+grid_Nghost_zmin(theGrid)+grid_Nghost_zmax(theGrid);
  int NUM_Q = grid_NUM_Q(theGrid);

  int i,j,k;
  for( k=0 ; k<N_z_withghost ; ++k ){
    for( i=0 ; i<N_r_withghost ; ++i ){
      for( j=0 ; j<grid_N_p(theGrid,i) ; ++j ){
        double phi = theCells[k][i][j].tiph;
        while( phi > 2.*M_PI ) phi -= 2.*M_PI;
        while( phi < 0.0 ) phi += 2.*M_PI;
        theCells[k][i][j].tiph = phi;
      }
    }
  }
}

void cell_copy(struct Cell ***theCells,struct Grid * theGrid){
  int N_r_withghost = grid_N_r(theGrid)+grid_Nghost_rmin(theGrid)+grid_Nghost_rmax(theGrid);
  int N_z_withghost = grid_N_z(theGrid)+grid_Nghost_zmin(theGrid)+grid_Nghost_zmax(theGrid);
  int NUM_Q = grid_NUM_Q(theGrid);
  int i,j,k,q;
  for( k=0 ; k<N_z_withghost ; ++k ){
    for( i=0 ; i<N_r_withghost ; ++i ){
      for( j=0 ; j<grid_N_p(theGrid,i) ; ++j ){
        for (q=0;q<NUM_Q;++q){
          theCells[k][i][j].RKcons[q] = theCells[k][i][j].cons[q];
        }
        theCells[k][i][j].RKtiph = theCells[k][i][j].tiph;
      }
    }
  }
}

void cell_adjust_RK_cons( struct Cell *** theCells , struct Grid * theGrid, double RK ){
  int N_r_withghost = grid_N_r(theGrid)+grid_Nghost_rmin(theGrid)+grid_Nghost_rmax(theGrid);
  int N_z_withghost = grid_N_z(theGrid)+grid_Nghost_zmin(theGrid)+grid_Nghost_zmax(theGrid);
  int NUM_Q = grid_NUM_Q(theGrid);

  int i,j,k,q;
  for( k=0 ; k<N_z_withghost ; ++k ){
    for( i=0 ; i<N_r_withghost ; ++i ){
      for( j=0 ; j<grid_N_p(theGrid,i) ; ++j ){
        struct Cell * c = &(theCells[k][i][j]);
        for( q=0 ; q<NUM_Q ; ++q ){
          c->cons[q] = (1.0-RK)*c->cons[q] + RK*c->RKcons[q];
        }
      }
    }
  }
}

void cell_update_phi( struct Cell *** theCells , struct Grid * theGrid, double RK , double dt ){
  int N_r_withghost = grid_N_r(theGrid)+grid_Nghost_rmin(theGrid)+grid_Nghost_rmax(theGrid);
  int N_z_withghost = grid_N_z(theGrid)+grid_Nghost_zmin(theGrid)+grid_Nghost_zmax(theGrid);

  int i,j,k;
  for( k=0 ; k<N_z_withghost ; ++k ){
    for( i=0 ; i<N_r_withghost ; ++i ){
      double rm = grid_r_faces(theGrid,i-1);
      double rp = grid_r_faces(theGrid,i);
      double r = .5*(rm+rp);
      for( j=0 ; j<grid_N_p(theGrid,i) ; ++j ){
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

void cell_update_dphi( struct Cell *** theCells,struct Grid * theGrid ){
  int N_r_withghost = grid_N_r(theGrid)+grid_Nghost_rmin(theGrid)+grid_Nghost_rmax(theGrid);
  int N_z_withghost = grid_N_z(theGrid)+grid_Nghost_zmin(theGrid)+grid_Nghost_zmax(theGrid);

  int i,j,k;
  for( k=0 ; k<N_z_withghost ; ++k ){
    for( i=0 ; i<N_r_withghost ; ++i ){
      for( j=0 ; j<grid_N_p(theGrid,i) ; ++j ){
        int jm = j-1;
        if( jm==-1 ) jm = grid_N_p(theGrid,i)-1;
        double dphi = theCells[k][i][j].tiph - theCells[k][i][jm].tiph;
        while( dphi > 2.*M_PI ) dphi -= 2.*M_PI;
        while( dphi < 0.0 ) dphi += 2.*M_PI;
        theCells[k][i][j].dphi = dphi;
      }
    }
  }
}



