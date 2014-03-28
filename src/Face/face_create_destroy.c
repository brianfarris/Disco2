#define FACE_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Sim.h"
#include "../Headers/Face.h"
#include "../Headers/Cell.h"
#include "../Headers/TimeStep.h"
#include "../Headers/header.h"

void addFace( struct Face * theFaces , int n , struct Cell * cL , struct Cell * cR , double r , double deltaL , double deltaR , double dphi , double tp , double deltaPerp ){
  theFaces[n].L = cL;
  theFaces[n].R = cR;
  theFaces[n].r = r;
  theFaces[n].deltaL = deltaL;
  theFaces[n].deltaR = deltaR;
  theFaces[n].dphi= dphi;
  theFaces[n].dA  = r*dphi*deltaPerp;
  theFaces[n].cm  = tp - .5*dphi;
} 

void build_jloop(int *pn,int i, int k,int rDir,int zDir,struct Cell *** theCells,struct Face * theFaces,struct Sim * theSim, int mode){

  double deltaL,deltaR,deltaPerp;
  double r;
  //printf("b1b1\n");
  if (rDir){
    deltaL = .5*(sim_FacePos(theSim,i,R_DIR)-sim_FacePos(theSim,i-1,R_DIR));
    deltaR = .5*(sim_FacePos(theSim,i+1,R_DIR)-sim_FacePos(theSim,i,R_DIR));
    r = sim_FacePos(theSim,i,R_DIR);

    double zp = sim_FacePos(theSim,k,Z_DIR);
    double zm = sim_FacePos(theSim,k-1,Z_DIR);
    deltaPerp = zp-zm;
  }
  //printf("b1b2\n");
  if (zDir){
    deltaL = .5*(sim_FacePos(theSim,k,Z_DIR)-sim_FacePos(theSim,k-1,Z_DIR));
    deltaR = .5*(sim_FacePos(theSim,k+1,Z_DIR)-sim_FacePos(theSim,k,Z_DIR));

    double rp=sim_FacePos(theSim,i,R_DIR);
    double rm = sim_FacePos(theSim,i-1,R_DIR);
    deltaPerp = rp-rm;
    r = .5*(rp+rm);
  }
  //printf("b1b3\n");

  int j,jp;
  double p0 = cell_tiph(cell_single(theCells,i,sim_N_p(theSim,i)-1,k));
  int jpmin=0;
  double dpmin=PHIMAX;
  //printf("b1b4\n");
  for( jp=0 ; jp<sim_N_p(theSim,i+rDir) ; ++jp ){
    double dp = cell_tiph(cell_single(theCells,i+rDir,jp,k+zDir))-p0;
    while( dp > PHIMAX ) dp -= PHIMAX;
    while( dp < 0.0 ) dp += PHIMAX;
    if( dpmin > dp ){
      dpmin = dp;
      jpmin = jp;
    }
  }
  //printf("b1b5\n");
  jp = jpmin;
  for( j=0 ; j<sim_N_p(theSim,i) ; ++j ){
    //printf("b1b6\n");
    int jm = j-1;
    if(jm<0) jm = sim_N_p(theSim,i)-1;
    double dp  = cell_tiph(cell_single(theCells,i+rDir,jp,k+zDir))-cell_tiph(cell_single(theCells,i,jm,k));
    while( dp > PHIMAX ) dp -= PHIMAX;
    while( dp < 0.0 ) dp += PHIMAX;
    //First figure out if cell+ covers all of cell-, 
    //if so create one face out of cell-.
    //printf("b1b7\n");
    if( cell_dphi(cell_single(theCells,i,j,k)) < dp ){
      //printf("b1b8\n");
      if ( mode==1 ){
        addFace(theFaces,*pn,cell_single(theCells,i,j,k),cell_single(theCells,i+rDir,jp,k+zDir),
            r,deltaL,deltaR,cell_dphi(cell_single(theCells,i,j,k)),cell_tiph(cell_single(theCells,i,j,k)),deltaPerp);
      }
      ++(*pn);
    }else{
      //printf("b1b9\n");
      //Otherwise, three steps:
      //Step A: face formed out of beginning of cell- and end of cell+. ++jp;
      if ( mode==1 ){
        addFace(theFaces,*pn,cell_single(theCells,i,j,k),cell_single(theCells,i+rDir,jp,k+zDir),
            r,deltaL,deltaR,dp,cell_tiph(cell_single(theCells,i+rDir,jp,k+zDir)),deltaPerp);
      }
      //printf("b1b10\n");
      ++(*pn);
      ++jp;
      if( jp == sim_N_p(theSim,i+rDir) ) {
        jp = 0;
      }
      //printf("b1b11\n");
      dp  = cell_tiph(cell_single(theCells,i,j,k))-cell_tiph(cell_single(theCells,i+rDir,jp,k+zDir));
      while( dp > PHIMAX/2. ) dp -= PHIMAX;
      //printf("b1b12\n");
      while( dp < -PHIMAX/2. ) dp += PHIMAX;
      //printf("b1b13\n");
      while( dp > 0.0 ){
        //Step B: (optional) all faces formed out of part of cell- and all of cell+. ++jp;
        if ( mode==1 ){
          addFace(theFaces,*pn,cell_single(theCells,i,j,k),cell_single(theCells,i+rDir,jp,k+zDir),
              r,deltaL,deltaR,cell_dphi(cell_single(theCells,i+rDir,jp,k+zDir)),cell_tiph(cell_single(theCells,i+rDir,jp,k+zDir)),deltaPerp);
        }
        //printf("b1b14\n");
        ++(*pn);
        ++jp;
        if( jp == sim_N_p(theSim,i+rDir) ) jp = 0;
        dp = cell_tiph(cell_single(theCells,i,j,k))-cell_tiph(cell_single(theCells,i+rDir,jp,k+zDir));
        //printf("b1b15\n");
        while( dp > PHIMAX/2. ) dp -= PHIMAX;
        //printf("b1b16\n");
        while( dp < -PHIMAX/2. ) dp += PHIMAX;
        //printf("b1b17\n");
      }
      //printf("b1b18\n");
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

struct Face *face_create(struct Cell *** theCells ,struct Sim *theSim, struct TimeStep * theTimeStep,int direction){

  struct Face * theFaces;

  if (direction==0){ // r-direction

    //printf("b1\n");
    //Count them first.  build_jloop with argument zero says "just count", doesn't create any faces.
    int i,k; 
    int n=0;
    for( i=0 ; i<sim_N(theSim,R_DIR)-1 ; ++i ){
      //printf("b1a\n");
      timestep_set_nri(theTimeStep,i,n);
      //printf("b1b\n");
      for( k=0 ; k<sim_N(theSim,Z_DIR) ; ++k ){
        build_jloop(&n,i,k,1,0,theCells,theFaces,theSim,0);
      }
      //printf("b1c\n");
    }
    //printf("b2\n");
    timestep_set_nri(theTimeStep,sim_N(theSim,R_DIR)-1,n);
    //timestep_set_Nfr(theTimeStep,theSim);
    //printf("b3\n");
    //allocate memory for array of Faces
    theFaces = (struct Face *) malloc( timestep_n(theTimeStep,sim_N(theSim,R_DIR)-1,R_DIR)*sizeof(struct Face) );
    //printf("b4\n");
    //now actually build the faces
    n=0;
    for( i=0 ; i<sim_N(theSim,R_DIR)-1 ; ++i ){
      for( k=0 ; k<sim_N(theSim,Z_DIR) ; ++k ){
        build_jloop(&n,i,k,1,0,theCells,theFaces,theSim,1);
      }
    }
    //printf("b5\n");
  }
  if (direction==1){ // z-direction

    if (sim_N_global(theSim,Z_DIR)>1){
      //Count them first.  build_jloop with argument zero says "just count", doesn't create any faces.
      int i,k;
      int n=0;
      for( k=0 ; k<sim_N(theSim,Z_DIR)-1 ; ++k ){
        timestep_set_nzk(theTimeStep,k,n);
        for( i=0 ; i<sim_N(theSim,R_DIR) ; ++i ){
          build_jloop(&n,i,k,0,1,theCells,theFaces,theSim,0);
        }
      }
      timestep_set_nzk(theTimeStep,sim_N(theSim,Z_DIR)-1,n);
      //timestep_set_Nfz(theTimeStep,theSim);

      //allocate memory for array of Faces
      theFaces = (struct Face *) malloc(timestep_n(theTimeStep,sim_N(theSim,Z_DIR)-1,Z_DIR)*sizeof(struct Face) );

      //now actually build the faces
      n=0;
      for( k=0 ; k<sim_N(theSim,Z_DIR)-1 ; ++k ){
        for( i=0 ; i<sim_N(theSim,R_DIR) ; ++i ){
          build_jloop(&n,i,k,0,1,theCells,theFaces,theSim,1);
        }
      }
    }
  }
  return(theFaces);
}

void face_destroy(struct Face * theFaces){
  free(theFaces);
}

