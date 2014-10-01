#define RIEMANN_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Riemann.h"
#include "../Headers/Sim.h"
#include "../Headers/Cell.h"
#include "../Headers/Face.h"
#include "../Headers/header.h"

void riemann_setup_rz(struct Riemann * theRiemann,struct Face * theFaces,struct Sim * theSim,int FaceNumber,int direction){
    theRiemann->n[direction]=1; // set
    int NUM_Q = sim_NUM_Q(theSim);
    double deltaL = face_deltaL(theFaces,FaceNumber);
    double deltaR = face_deltaR(theFaces,FaceNumber);
    theRiemann->cL = face_L_pointer(theFaces,FaceNumber);
    theRiemann->cR = face_R_pointer(theFaces,FaceNumber);
    double pL = cell_tiph(theRiemann->cL) - .5*cell_dphi(theRiemann->cL);
    double pR = cell_tiph(theRiemann->cR) - .5*cell_dphi(theRiemann->cR);   
    double dpL =  face_cm(theFaces,FaceNumber) - pL;
    double dpR = -face_cm(theFaces,FaceNumber) + pR;
    while( dpL >  PHIMAX/2. ) dpL -= PHIMAX;
    while( dpL < -PHIMAX/2. ) dpL += PHIMAX;
    while( dpR >  PHIMAX/2. ) dpR -= PHIMAX;
    while( dpR < -PHIMAX/2. ) dpR += PHIMAX;
    dpL = dpL;
    dpR = dpR;
    theRiemann->r = face_r(theFaces,FaceNumber);
    theRiemann->dA = face_dA(theFaces,FaceNumber);
    theRiemann->cm = face_cm(theFaces,FaceNumber);
    if (direction==0){
        theRiemann->r_cell_L = theRiemann->r-deltaL;
        theRiemann->r_cell_R = theRiemann->r+deltaR;
    } else{
        theRiemann->r_cell_L = theRiemann->r;
        theRiemann->r_cell_R = theRiemann->r;
    }

    int q;
    for (q=0;q<NUM_Q;++q){
        theRiemann->primL[q] = cell_prim(theRiemann->cL,q) + cell_grad(theRiemann->cL,q)*deltaL + cell_gradp(theRiemann->cL,q)*dpL;
        theRiemann->primR[q] = cell_prim(theRiemann->cR,q) - cell_grad(theRiemann->cR,q)*deltaR - cell_gradp(theRiemann->cR,q)*dpR;
    }
}

void riemann_setup_p(struct Riemann * theRiemann,struct Cell *** theCells,struct Sim * theSim,int i,int j_low,int k,int direction){
    theRiemann->n[direction]=1; // set
    int NUM_Q = sim_NUM_Q(theSim);

    int j_hi;
    if (j_low == sim_N_p(theSim,i)-1){
        j_hi = 0;
    } else{
        j_hi = j_low+1;
    }
    theRiemann->cL = cell_single(theCells,i,j_low,k);
    theRiemann->cR = cell_single(theCells,i,j_hi ,k);
    double dpL = cell_dphi(theRiemann->cL);
    double dpR = cell_dphi(theRiemann->cR);
    double zm = sim_FacePos(theSim,k-1,Z_DIR);
    double zp = sim_FacePos(theSim,k,Z_DIR);
    double dz = zp-zm;
    double rm = sim_FacePos(theSim,i-1,R_DIR);
    double rp = sim_FacePos(theSim,i,R_DIR);
    double dr = rp-rm;
    double r = .5*(rp+rm);
    theRiemann->dA = dr*dz;
    theRiemann->r = r; 
    theRiemann->cm = cell_tiph(theRiemann->cL);

    theRiemann->r_cell_L = r;
    theRiemann->r_cell_R = r;

    int q;
    for (q=0;q<NUM_Q;++q){
        theRiemann->primL[q] = cell_prim(theRiemann->cL,q) + 0.5*cell_gradp(theRiemann->cL,q)*dpL;
        theRiemann->primR[q] = cell_prim(theRiemann->cR,q) - 0.5*cell_gradp(theRiemann->cR,q)*dpR;
    }

}

