#define FACE_PRIVATE_DEFS
#include <stdlib.h>
#include <stdIO.h>
#include <math.h>
#include "../Headers/Grid.h"
#include "../Headers/Face.h"
#include "../Headers/Cell.h"
#include "../Headers/header.h"

void addFace( struct Face * theFaces , int n , struct Cell * cL , struct Cell * cR , double r , double deltaL , double deltaR , double dphi , double tp , double dz ){
  theFaces[n].L = cL;
  theFaces[n].R = cR;
  theFaces[n].r = r;
  theFaces[n].deltaL = deltaL;
  theFaces[n].deltaR = deltaR;
  theFaces[n].dphi= dphi;
  theFaces[n].dA  = r*dphi*dz;
  theFaces[n].cm  = tp - .5*dphi;
} 

void face_build_r( struct Cell *** theCells , struct Face * theFaces , int * nri , int mode , struct Grid *theGrid){
  int N_r_withghost = grid_N_r(theGrid)+grid_Nghost_rmin(theGrid)+grid_Nghost_rmax(theGrid);
  int N_z_withghost = grid_N_z(theGrid)+grid_Nghost_zmin(theGrid)+grid_Nghost_zmax(theGrid);

  int i,j,k; 
  int n=0;
  for( i=0 ; i<N_r_withghost-1 ; ++i ){
    if( mode == 0 ) nri[i] = n;
    for( k=0 ; k<N_z_withghost ; ++k ){
      double deltaL,deltaR;
      deltaL = .5*(grid_r_faces(theGrid,i)-grid_r_faces(theGrid,i-1));
      deltaR = .5*(grid_r_faces(theGrid,i+1)-grid_r_faces(theGrid,i));

      double zp = grid_z_faces(theGrid,k);
      double zm = grid_z_faces(theGrid,k-1);
      double dz = zp-zm;

      int jp;
      double p0 = cell_tiph(theCells,i,grid_N_p(theGrid,i)-1,k);
      int jpmin=0;
      double dpmin=2.*M_PI;
      for( jp=0 ; jp<grid_N_p(theGrid,i+1) ; ++jp ){
        double dp = cell_tiph(theCells,i+1,jp,k) - p0;
        while( dp > 2.*M_PI ) dp -= 2.*M_PI;
        while( dp < 0.0 ) dp += 2.*M_PI;
        if( dpmin > dp ){
          dpmin = dp;
          jpmin = jp;
        }
      }
      jp = jpmin;
      for( j=0 ; j<grid_N_p(theGrid,i) ; ++j ){
        int jm = j-1;
        if(jm<0) jm = grid_N_p(theGrid,i)-1;
        double dp  = cell_tiph(theCells,i+1,jp,k)-cell_tiph(theCells,i,jm,k);
        while( dp > 2.*M_PI ) dp -= 2.*M_PI;
        while( dp < 0.0 ) dp += 2.*M_PI;
        //First figure out if cell+ covers all of cell-, 
        //if so create one face out of cell-.
        if( cell_dphi(theCells,i,j,k) < dp ){
          if( mode==1 ){
            addFace(theFaces,n,cell_pointer(theCells,i,j,k),cell_pointer(theCells,i+1,jp,k),
                grid_r_faces(theGrid,i),deltaL,deltaR,cell_dphi(theCells,i,j,k),cell_tiph(theCells,i,j,k),dz);
          }
          ++n;
        }else{
          //Otherwise, three steps:
          //Step A: face formed out of beginning of cell- and end of cell+. ++jp;
          if( mode==1 ){
            addFace(theFaces,n,cell_pointer(theCells,i,j,k),cell_pointer(theCells,i+1,jp,k),
                grid_r_faces(theGrid,i),deltaL,deltaR,dp,cell_tiph(theCells,i+1,jp,k),dz);
          }
          ++n;
          ++jp;
          if( jp == grid_N_p(theGrid,i+1) ) jp = 0;
          dp  = cell_tiph(theCells,i,j,k) - cell_tiph(theCells,i+1,jp,k);
          while( dp > M_PI ) dp -= 2.*M_PI;
          while( dp < -M_PI ) dp += 2.*M_PI;
          while( dp > 0.0 ){
            //Step B: (optional) all faces formed out of part of cell- and all of cell+. ++jp;
            if( mode==1 ){
              addFace(theFaces,n,cell_pointer(theCells,i,j,k),cell_pointer(theCells,i+1,jp,k),
                  grid_r_faces(theGrid,i),deltaL,deltaR,cell_dphi(theCells,i+1,jp,k),cell_tiph(theCells,i+1,jp,k) ,dz);
            }
            ++n;
            ++jp;
            if( jp == grid_N_p(theGrid,i+1) ) jp = 0;
            dp = cell_tiph(theCells,i,j,k)-cell_tiph(theCells,i+1,jp,k);
            while( dp > M_PI ) dp -= 2.*M_PI;
            while( dp < -M_PI ) dp += 2.*M_PI;
          }
          dp = cell_dphi(theCells,i+1,jp,k)+dp;
          //Step C: face formed out of end of cell- and beginning of cell+.
          if( mode==1 ){
            addFace(theFaces,n,cell_pointer(theCells,i,j,k),cell_pointer(theCells,i+1,jp,k),
                grid_r_faces(theGrid,i),deltaL,deltaR,dp,cell_tiph(theCells,i,j,k),dz);
          }
          ++n;
        }
      }
    }
  }
  if( mode==0 ) nri[N_r_withghost-1] = n;
 
}

void face_build_z( struct Cell *** theCells , struct Face * theFaces, int * nzk , int mode , struct Grid *theGrid){
  int N_r_withghost = grid_N_r(theGrid)+grid_Nghost_rmin(theGrid)+grid_Nghost_rmax(theGrid);
  int N_z_withghost = grid_N_z(theGrid)+grid_Nghost_zmin(theGrid)+grid_Nghost_zmax(theGrid);

  int i,j,k;
  int n=0;
  for( k=0 ; k<N_z_withghost-1 ; ++k ){
    if( mode == 0 ) nzk[k] = n;
    for( i=0 ; i<N_r_withghost ; ++i ){

      double rp=grid_r_faces(theGrid,i);
      double rm = grid_r_faces(theGrid,i-1);
      double dr = rp-rm;
      double r = .5*(rp+rm);

      double dzL,dzR;
      dzL = .5*(grid_z_faces(theGrid,k)-grid_z_faces(theGrid,k-1));
      dzR = .5*(grid_z_faces(theGrid,k+1)-grid_z_faces(theGrid,k));

      int jp;
      double p0 = cell_tiph(theCells,i,grid_N_p(theGrid,i)-1,k);
      int jpmin=0;
      double dpmin=2.*M_PI;
      for( jp=0 ; jp<grid_N_p(theGrid,i) ; ++jp ){
        double dp = cell_tiph(theCells,i,jp,k+1)-p0;
        while( dp > 2.*M_PI ) dp -= 2.*M_PI;
        while( dp < 0.0 ) dp += 2.*M_PI;
        if( dpmin > dp ){
          dpmin = dp;
          jpmin = jp;
        }
      }
      jp = jpmin;
      for( j=0 ; j<grid_N_p(theGrid,i) ; ++j ){
        int jm = j-1;
        if(jm<0) jm = grid_N_p(theGrid,i)-1;
        double dp  = cell_tiph(theCells,i,jp,k+1)-cell_tiph(theCells,i,jm,k);
        while( dp > 2.*M_PI ) dp -= 2.*M_PI;
        while( dp < 0.0 ) dp += 2.*M_PI;
        //First figure out if cell+ covers all of cell-, 
        //if so create one face out of cell-.
        if( cell_dphi(theCells,i,j,k) < dp ){
          if(mode==1){
            addFace(theFaces,n,cell_pointer(theCells,i,j,k),cell_pointer(theCells,i,jp,k+1),
                r,dzL,dzR,cell_dphi(theCells,i,j,k),cell_tiph(theCells,i,j,k),dr);
          }
          ++n;
        }else{
          //Otherwise, three steps:
          //Step A: face formed out of beginning of cell- and end of cell+. ++jp;
          if (mode==1 ){
            addFace(theFaces,n,cell_pointer(theCells,i,j,k),cell_pointer(theCells,i,jp,k+1),
                r,dzL,dzR,dp,cell_tiph(theCells,i,jp,k+1),dr);
          }
          ++n;
          ++jp;
          if( jp == grid_N_p(theGrid,i) ) jp = 0;
          dp  = cell_tiph(theCells,i,j,k) - cell_tiph(theCells,i,jp,k+1);
          while( dp > M_PI ) dp -= 2.*M_PI;
          while( dp < -M_PI ) dp += 2.*M_PI;
          while( dp > 0.0 ){
            //Step B: (optional) all faces formed out of part of cell- and all of cell+. ++jp;
            if(mode==1){
              addFace(theFaces,n,cell_pointer(theCells,i,j,k),cell_pointer(theCells,i,jp,k+1),
                  r,dzL,dzR,cell_dphi(theCells,i,jp,k+1),cell_tiph(theCells,i,jp,k+1),dr);
            }
            ++n;
            ++jp;
            if( jp == grid_N_p(theGrid,i) ) jp = 0;
            dp  = cell_tiph(theCells,i,j,k) - cell_tiph(theCells,i,jp,k+1);
            while( dp > M_PI ) dp -= 2.*M_PI;
            while( dp < -M_PI ) dp += 2.*M_PI;
          }
          dp = cell_dphi(theCells,i,jp,k+1)+dp;
          //Step C: face formed out of end of cell- and beginning of cell+.
          if( mode==1 ){
            addFace(theFaces,n,cell_pointer(theCells,i,j,k),cell_pointer(theCells,i,jp,k+1),
                r,dzL,dzR,dp,cell_tiph(theCells,i,j,k),dr);
          }
          ++n;
        }
      }
    }
  }
  if( mode==0 ) nzk[N_z_withghost-1] = n;
}

struct Face *face_create_r(struct Cell *** theCells ,struct Grid *theGrid, int *pNfr, int *nri) {
  int N_r_withghost = grid_N_r(theGrid)+grid_Nghost_rmin(theGrid)+grid_Nghost_rmax(theGrid);
  struct Face * theFaces_r;
  //int nri[N_r_withghost];
  //Count them first.  Buildfaces with argument zero says "just count", doesn't create any faces.
  face_build_r( theCells , NULL , nri , 0 ,theGrid);
   *pNfr = nri[N_r_withghost-1];
  int Nfr = *pNfr;
  theFaces_r = (struct Face *) malloc( Nfr*sizeof(struct Face) );
  face_build_r( theCells , theFaces_r , nri , 1 ,theGrid);
   return(theFaces_r);
}

struct Face *face_create_z(struct Cell *** theCells ,struct Grid *theGrid, int *pNfz, int *nzk) {
  int N_z_withghost = grid_N_z(theGrid)+grid_Nghost_zmin(theGrid)+grid_Nghost_zmax(theGrid);

  struct Face * theFaces_z;
  //int nzk[N_z_withghost];
  //Count them first.  Buildfaces with argument zero says "just count", doesn't create any faces.
  face_build_z( theCells , NULL , nzk , 0 ,theGrid);
  *pNfz = nzk[N_z_withghost-1];
  int Nfz = *pNfz;
  theFaces_z = (struct Face *) malloc( Nfz*sizeof(struct Face) );
  face_build_z( theCells , theFaces_z , nzk , 1 ,theGrid);
  return(theFaces_z);
}

void face_destroy(struct Face * theFaces){
  free(theFaces);
}

