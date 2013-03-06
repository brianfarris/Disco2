#define FACE_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Grid.h"
#include "../Headers/Face.h"
#include "../Headers/Cell.h"
#include "../Headers/TimeStep.h"
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

void build_jloop(int *pn,int i, int k,int rDir,int zDir,struct Cell *** theCells,struct Face * theFaces,struct Grid * theGrid, int mode){

  double deltaL,deltaR,deltaPerp;
  double r;
  if (rDir){
    deltaL = .5*(grid_r_faces(theGrid,i)-grid_r_faces(theGrid,i-1));
    deltaR = .5*(grid_r_faces(theGrid,i+1)-grid_r_faces(theGrid,i));
    r = grid_r_faces(theGrid,i);

    double zp = grid_z_faces(theGrid,k);
    double zm = grid_z_faces(theGrid,k-1);
    deltaPerp = zp-zm;
  }
  if (zDir){
    deltaL = .5*(grid_z_faces(theGrid,k)-grid_z_faces(theGrid,k-1));
    deltaR = .5*(grid_z_faces(theGrid,k+1)-grid_z_faces(theGrid,k));

    double rp=grid_r_faces(theGrid,i);
    double rm = grid_r_faces(theGrid,i-1);
    deltaPerp = rp-rm;
    r = .5*(rp+rm);
  }

  int j,jp;
  double p0 = cell_tiph(cell_single(theCells,i,grid_N_p(theGrid,i)-1,k));
  int jpmin=0;
  double dpmin=2.*M_PI;
  for( jp=0 ; jp<grid_N_p(theGrid,i+rDir) ; ++jp ){
    double dp = cell_tiph(cell_single(theCells,i+rDir,jp,k+zDir))-p0;
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
    double dp  = cell_tiph(cell_single(theCells,i+rDir,jp,k+zDir))-cell_tiph(cell_single(theCells,i,jm,k));
    while( dp > 2.*M_PI ) dp -= 2.*M_PI;
    while( dp < 0.0 ) dp += 2.*M_PI;
    //First figure out if cell+ covers all of cell-, 
    //if so create one face out of cell-.
    if( cell_dphi(cell_single(theCells,i,j,k)) < dp ){
      if ( mode==1 ){
        addFace(theFaces,*pn,cell_single(theCells,i,j,k),cell_single(theCells,i+rDir,jp,k+zDir),
            r,deltaL,deltaR,cell_dphi(cell_single(theCells,i,j,k)),cell_tiph(cell_single(theCells,i,j,k)),deltaPerp);
      }
      ++(*pn);
    }else{
      //Otherwise, three steps:
      //Step A: face formed out of beginning of cell- and end of cell+. ++jp;
      if ( mode==1 ){
        addFace(theFaces,*pn,cell_single(theCells,i,j,k),cell_single(theCells,i+rDir,jp,k+zDir),
            r,deltaL,deltaR,dp,cell_tiph(cell_single(theCells,i+rDir,jp,k+zDir)),deltaPerp);
      }
      ++(*pn);
      ++jp;
      if( jp == grid_N_p(theGrid,i+rDir) ) jp = 0;
      dp  = cell_tiph(cell_single(theCells,i,j,k))-cell_tiph(cell_single(theCells,i+rDir,jp,k+zDir));
      while( dp > M_PI ) dp -= 2.*M_PI;
      while( dp < -M_PI ) dp += 2.*M_PI;
      while( dp > 0.0 ){
        //Step B: (optional) all faces formed out of part of cell- and all of cell+. ++jp;
        if ( mode==1 ){
          addFace(theFaces,*pn,cell_single(theCells,i,j,k),cell_single(theCells,i+rDir,jp,k+zDir),
              r,deltaL,deltaR,cell_dphi(cell_single(theCells,i+rDir,jp,k+zDir)),cell_tiph(cell_single(theCells,i+rDir,jp,k+zDir)),deltaPerp);
        }
        ++(*pn);
        ++jp;
        if( jp == grid_N_p(theGrid,i+rDir) ) jp = 0;
        dp = cell_tiph(cell_single(theCells,i,j,k))-cell_tiph(cell_single(theCells,i+rDir,jp,k+zDir));
        while( dp > M_PI ) dp -= 2.*M_PI;
        while( dp < -M_PI ) dp += 2.*M_PI;
      }
      dp = cell_dphi(cell_single(theCells,i+rDir,jp,k+zDir))+dp;
      //Step C: face formed out of end of cell- and beginning of cell+.
      if ( mode==1 ){
        addFace(theFaces,*pn,cell_single(theCells,i,j,k),cell_single(theCells,i+rDir,jp,k+zDir),
            r,deltaL,deltaR,dp,cell_tiph(cell_single(theCells,i,j,k)),deltaPerp);
      }
      ++(*pn);
    }
  }
}

//try to consolidate using function pointers
struct Face *face_create(struct Cell *** theCells ,struct Grid *theGrid, struct TimeStep * theTimeStep,int direction){
  int N_r_withghost = grid_N_r(theGrid)+grid_Nghost_rmin(theGrid)+grid_Nghost_rmax(theGrid);
  int N_z_withghost = grid_N_z(theGrid)+grid_Nghost_zmin(theGrid)+grid_Nghost_zmax(theGrid);

  struct Face * theFaces;

  if (direction==0){ // r-direction

    //Count them first.  build_jloop with argument zero says "just count", doesn't create any faces.
    int i,k; 
    int n=0;
    for( i=0 ; i<N_r_withghost-1 ; ++i ){
      timestep_nri(theTimeStep)[i] = n;
      for( k=0 ; k<N_z_withghost ; ++k ){
        build_jloop(&n,i,k,1,0,theCells,theFaces,theGrid,0);
      }
    }
    timestep_nri(theTimeStep)[N_r_withghost-1] = n; 
    timestep_set_Nfr(theTimeStep,theGrid);

    //allocate memory for array of Faces
    theFaces = (struct Face *) malloc( timestep_Nfr(theTimeStep)*sizeof(struct Face) );

    //now actually build the faces
    n=0;
    for( i=0 ; i<N_r_withghost-1 ; ++i ){
      for( k=0 ; k<N_z_withghost ; ++k ){
        build_jloop(&n,i,k,1,0,theCells,theFaces,theGrid,1);
      }
    }
  }

  if (direction==1){ // z-direction

    if (grid_N_z_global(theGrid)>1){
      //Count them first.  build_jloop with argument zero says "just count", doesn't create any faces.
      int i,k;
      int n=0;
      for( k=0 ; k<N_z_withghost-1 ; ++k ){
        timestep_nzk(theTimeStep)[k] = n;
        for( i=0 ; i<N_r_withghost ; ++i ){
          build_jloop(&n,i,k,0,1,theCells,theFaces,theGrid,0);
        }
      }
      timestep_nzk(theTimeStep)[N_z_withghost-1] = n;
      timestep_set_Nfz(theTimeStep,theGrid);

      //allocate memory for array of Faces
      theFaces = (struct Face *) malloc( timestep_Nfz(theTimeStep)*sizeof(struct Face) );

      //now actually build the faces
      n=0;
      for( k=0 ; k<N_z_withghost-1 ; ++k ){
        for( i=0 ; i<N_r_withghost ; ++i ){
          build_jloop(&n,i,k,0,1,theCells,theFaces,theGrid,1);
        }
      }
    }
  }
  return(theFaces);
}

void face_destroy(struct Face * theFaces){
  free(theFaces);
}

