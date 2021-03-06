#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Sim.h"
#include "../Headers/Face.h"
#include "../Headers/GravMass.h"
#include "../Headers/TimeStep.h"
#include "../Headers/MPIsetup.h"
#include "../Headers/header.h"

void cell_plm_rz( struct Cell *** theCells ,struct Sim *theSim, struct Face * theFaces , struct TimeStep * theTimeStep , struct MPIsetup * theMPIsetup ,  int direction ){
  int NUM_Q = sim_NUM_Q(theSim);
  double PLM = sim_PLM(theSim);

  int Nf = timestep_n(theTimeStep,sim_N(theSim,direction)-1,direction); // number of faces
  
  int i,j,k,q;
  for( k=0 ; k<sim_N(theSim,Z_DIR) ; ++k ){
    for( i=0 ; i<sim_N(theSim,R_DIR) ; ++i ){
      for( j=0 ; j<sim_N_p(theSim,i) ; ++j ){
        for( q=0 ; q<NUM_Q ; ++q ){
          theCells[k][i][j].grad[q] = 0.0;
        }
      }
    }
  }

  int n;
  for( n=0 ; n<Nf ; ++n ){
    struct Cell * cL = face_L_pointer(theFaces,n);
    struct Cell * cR = face_R_pointer(theFaces,n);
    double deltaL = face_deltaL(theFaces,n);
    double deltaR = face_deltaR(theFaces,n);
    double pL = cL->tiph - .5*cL->dphi;
    double pR = cR->tiph - .5*cR->dphi;
    double dpL = face_cm(theFaces,n) - pL;
    while( dpL >  PHIMAX/2. ) dpL -= PHIMAX;
    while( dpL < -PHIMAX/2. ) dpL += PHIMAX;
    double dpR = pR - face_cm(theFaces,n);
    while( dpR >  PHIMAX/2. ) dpR -= PHIMAX;
    while( dpR < -PHIMAX/2. ) dpR += PHIMAX;
    double dA  = face_dA(theFaces,n);
    double face_rad = face_r(theFaces,n);
    for( q=0 ; q<NUM_Q ; ++q ){
      double WL = cL->prim[q] + dpL*cL->gradp[q];
      double WR = cR->prim[q] - dpR*cR->gradp[q];
      double S = (WR-WL)/(deltaR+deltaL);
      cL->grad[q] += S*dA;
      cR->grad[q] += S*dA; 
    }
  }
  for( k=0 ; k<sim_N(theSim,Z_DIR) ; ++k ){
    double zm = sim_FacePos(theSim,k-1,Z_DIR);
    double zp = sim_FacePos(theSim,k,Z_DIR);
    double dz = zp-zm;
    for( i=0 ; i<sim_N(theSim,R_DIR) ; ++i ){
      double rm = sim_FacePos(theSim,i-1,R_DIR);
      double rp = sim_FacePos(theSim,i,R_DIR);
      for( j=0 ; j<sim_N_p(theSim,i) ; ++j ){
        double dp = theCells[k][i][j].dphi;
        double dAtot;
        if( direction==R_DIR ){
          if (mpisetup_check_rin_bndry(theMPIsetup) && i==0){
            dAtot = dz*rp*dp;
          } else if(mpisetup_check_rout_bndry(theMPIsetup) && i==sim_N(theSim,R_DIR)-1){
            dAtot = dz*rm*dp; 
          }else {
            dAtot = dz*(rp+rm)*dp;
          }
        } else{
          dAtot = 2.*.5*(rp*rp-rm*rm)*dp;
        }
        for( q=0 ; q<NUM_Q ; ++q ){
          theCells[k][i][j].grad[q] /= dAtot;
        }
      }
    }
  }
  for( n=0 ; n<Nf ; ++n ){
    struct Cell * cL = face_L_pointer(theFaces,n);
    struct Cell * cR = face_R_pointer(theFaces,n);
    double deltaL = face_deltaL(theFaces,n);
    double deltaR = face_deltaR(theFaces,n);
    double pL = cL->tiph - .5*cL->dphi;
    double pR = cR->tiph - .5*cR->dphi;
    double dpL = face_cm(theFaces,n) - pL;
    while( dpL >  PHIMAX/2. ) dpL -= PHIMAX;
    while( dpL < -PHIMAX/2. ) dpL += PHIMAX;
    double dpR = pR - face_cm(theFaces,n);
    while( dpR >  PHIMAX/2. ) dpR -= PHIMAX;
    while( dpR < -PHIMAX/2. ) dpR += PHIMAX;

    for( q=0 ; q<NUM_Q ; ++q ){
      double WL = cL->prim[q] + dpL*cL->gradp[q];
      double WR = cR->prim[q] - dpR*cR->gradp[q];

      double S = (WR-WL)/(deltaR+deltaL);
      double SL = cL->grad[q];
      double SR = cR->grad[q];
      if( S*SL < 0.0 ){
        cL->grad[q] =0.0;
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

void cell_plm_p( struct Cell *** theCells ,struct Sim * theSim){
  int NUM_Q = sim_NUM_Q(theSim);
  double PLM = sim_PLM(theSim);

  int Qmin = 0;
  int Qmax = NUM_Q;

  int i,j,k,q;
  for( k=0 ; k<sim_N(theSim,Z_DIR) ; ++k ){
    for( i=0 ; i<sim_N(theSim,R_DIR) ; ++i ){
      for( j=0 ; j<sim_N_p(theSim,i) ; ++j ){
        struct Cell * c = &(theCells[k][i][j]);
        struct Cell * cL;
        struct Cell * cR;
        if( j==0 ){
          cL = &( theCells[k][i][sim_N_p(theSim,i)-1] );
          cR = &( theCells[k][i][j+1] );
        }else if(j==sim_N_p(theSim,i)-1){
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
        double rp = sim_FacePos(theSim,i,R_DIR);
        double rm = sim_FacePos(theSim,i-1,R_DIR);
        double r = .5*(rp+rm);
      }
    }
  }   
}


