#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Grid.h"
#include "../Headers/Face.h"
#include "../Headers/GravMass.h"
#include "../Headers/header.h"

void cell_plm_rz( struct Cell *** theCells ,struct Grid *theGrid, struct Face * theFaces , int Nf , int rz ){
  int N_r_withghost = grid_N_r(theGrid)+grid_Nghost_rmin(theGrid)+grid_Nghost_rmax(theGrid);
  int N_z_withghost = grid_N_z(theGrid)+grid_Nghost_zmin(theGrid)+grid_Nghost_zmax(theGrid);
  int NUM_Q = grid_NUM_Q(theGrid);
  double PLM = grid_PLM(theGrid);

  int i,j,k,q;

  for( k=0 ; k<N_z_withghost ; ++k ){
    for( i=0 ; i<N_r_withghost ; ++i ){
      for( j=0 ; j<grid_N_p(theGrid,i) ; ++j ){
        for( q=0 ; q<NUM_Q ; ++q ){
          theCells[k][i][j].grad[q] = 0.0;
        }
      }
    }
  }

  int n;
  for( n=0 ; n<Nf ; ++n ){
    struct Face * f  = face_pointer(theFaces,n);
    struct Cell * cL = face_L_pointer(theFaces,n);
    struct Cell * cR = face_R_pointer(theFaces,n);
    double deltaL = face_deltaL(f);
    double deltaR = face_deltaR(f);
    double pL = cL->tiph - .5*cL->dphi;
    double pR = cR->tiph - .5*cR->dphi;
    double dpL = face_cm(f) - pL;
    while( dpL >  M_PI ) dpL -= 2.*M_PI;
    while( dpL < -M_PI ) dpL += 2.*M_PI;
    double dpR = pR - face_cm(f);
    while( dpR >  M_PI ) dpR -= 2.*M_PI;
    while( dpR < -M_PI ) dpR += 2.*M_PI;
    double dA  = face_dA(f);
    for( q=0 ; q<NUM_Q ; ++q ){
      double WL = cL->prim[q] + dpL*cL->gradp[q];
      double WR = cR->prim[q] - dpR*cR->gradp[q];

      double S = (WR-WL)/(deltaR+deltaL);
      cL->grad[q] += S*dA;
      cR->grad[q] += S*dA; 
    }
  }
  for( k=0 ; k<N_z_withghost ; ++k ){
    double zm = grid_z_faces(theGrid,k-1);
    double zp = grid_z_faces(theGrid,k);
    double dz = zp-zm;
    for( i=0 ; i<N_r_withghost ; ++i ){
      double rm = grid_r_faces(theGrid,i-1);
      double rp = grid_r_faces(theGrid,i);
      for( j=0 ; j<grid_N_p(theGrid,i) ; ++j ){
        double dp = theCells[k][i][j].dphi;
        double dAtot;
        if( rz==0 ) dAtot = dz*(rp+rm)*dp; else dAtot = 2.*.5*(rp*rp-rm*rm)*dp;
        for( q=0 ; q<NUM_Q ; ++q ){
          theCells[k][i][j].grad[q] /= dAtot;
        }
      }
    }
  }

  for( n=0 ; n<Nf ; ++n ){
    struct Face * f  = face_pointer(theFaces,n);
    struct Cell * cL = face_L_pointer(theFaces,n);
    struct Cell * cR = face_R_pointer(theFaces,n);
    double deltaL = face_deltaL(f);
    double deltaR = face_deltaR(f);
    double pL = cL->tiph - .5*cL->dphi;
    double pR = cR->tiph - .5*cR->dphi;
    double dpL = face_cm(f) - pL;
    while( dpL >  M_PI ) dpL -= 2.*M_PI;
    while( dpL < -M_PI ) dpL += 2.*M_PI;
    double dpR = pR - face_cm(f);
    while( dpR >  M_PI ) dpR -= 2.*M_PI;
    while( dpR < -M_PI ) dpR += 2.*M_PI;
    for( q=0 ; q<NUM_Q ; ++q ){
      double WL = cL->prim[q] + dpL*cL->gradp[q];
      double WR = cR->prim[q] - dpR*cR->gradp[q];

      double S = (WR-WL)/(deltaR+deltaL);
      double SL = cL->grad[q];
      double SR = cR->grad[q];
      if( S*SL < 0.0 ){
        cL->grad[q] = 0.0;
      }else if( fabs(PLM*S) < fabs(SL) ){
        cL->grad[q] = PLM*S;
      }
      if( S*SR < 0.0 ){
        cR->grad[q] = 0.0;
      }else if( fabs(PLM*S) < fabs(SR) ){
        cR->grad[q] = PLM*S;
      }
    }
  }
}

void cell_plm_p( struct Cell *** theCells ,struct Grid * theGrid){
  int N_r_withghost = grid_N_r(theGrid)+grid_Nghost_rmin(theGrid)+grid_Nghost_rmax(theGrid);
  int N_z_withghost = grid_N_z(theGrid)+grid_Nghost_zmin(theGrid)+grid_Nghost_zmax(theGrid);
  int NUM_Q = grid_NUM_Q(theGrid);
  double PLM = grid_PLM(theGrid);

  int Qmin = 0;
  int Qmax = NUM_Q;

  int i,j,k,q;
  for( k=0 ; k<N_z_withghost ; ++k ){
    for( i=0 ; i<N_r_withghost ; ++i ){
      for( j=0 ; j<grid_N_p(theGrid,i) ; ++j ){
        struct Cell * c = &(theCells[k][i][j]);
        struct Cell * cL;
        struct Cell * cR;
        if( j==0 ){
          cL = &( theCells[k][i][grid_N_p(theGrid,i)-1] );
          cR = &( theCells[k][i][j+1] );
        }else if(j==grid_N_p(theGrid,i)-1){
          cL = &( theCells[k][i][j-1] );
          cR = &( theCells[k][i][0] );
        }else{
          cL = &( theCells[k][i][j-1] );
          cR = &( theCells[k][i][j+1] );
        }

        double dpL = .5*( c->dphi + cL->dphi );
        double dpR = .5*( c->dphi + cR->dphi );

        for( q=Qmin ; q<Qmax ; ++q ){
          double sL = ( c->prim[q]  - cL->prim[q] )/dpL;
          double sR = ( cR->prim[q] - c->prim[q]  )/dpR;
          double sM = ( cR->prim[q] - cL->prim[q] )/(dpL+dpR);

          double s = PLM*sL;
          if( fabs(s) > PLM*fabs(sR) ) s = PLM*sR;
          if( fabs(s) > fabs(sM) ) s = sM;
          if( sL*sR < 0.0 ) s = 0.0;
          c->gradp[q] = s;

        } 
      }
    }
  }   

}


