#define IO_PRIVATE_DEFS
#include <string.h>
#include <stdlib.h>
#include "../Headers/Sim.h"
#include "../Headers/Cell.h"
#include "../Headers/GravMass.h"
#include "hdf5.h"
#include "../Headers/IO.h"
#include "../Headers/header.h"

void io_allocbuf(struct IO *theIO,struct Sim *theSim){
  int NUM_Q = sim_NUM_Q(theSim);
  theIO->buffer = malloc(sizeof(double *)*sim_Ncells(theSim));
  theIO->buffer[0] = malloc(sizeof(double )*sim_Ncells(theSim)*(NUM_Q+3));
  int i;
  for (i=0;i<sim_Ncells(theSim);++i){
    theIO->buffer[i]= theIO->buffer[0]+i*(NUM_Q+3);
  }
}

void io_deallocbuf(struct IO *theIO){
  free(theIO->buffer[0]);
  free(theIO->buffer);
}

void io_setbuf(struct IO *theIO,struct Cell ***theCells,struct Sim *theSim,struct GravMass * theGravMasses){
  int NUM_Q = sim_NUM_Q(theSim);
  int i,j,k,q;
  int index=0;
  for (k=0; k<sim_N(theSim,Z_DIR); k++) {
    for (i=0; i<sim_N(theSim,R_DIR);i++){
      for (j=0; j<sim_N_p(theSim,i);j++){
        theIO->buffer[index][0] = cell_tiph(cell_single(theCells,i,j,k));
        theIO->buffer[index][1] = 0.5*(sim_FacePos(theSim,i-1,R_DIR)+sim_FacePos(theSim,i,R_DIR));
        theIO->buffer[index][2] = 0.5*(sim_FacePos(theSim,k-1,Z_DIR)+sim_FacePos(theSim,k,Z_DIR));
        for (q=0;q<NUM_Q;q++){
          theIO->buffer[index][q+3] = cell_prim(cell_single(theCells,i,j,k),q); 
        }
        index += 1;
      }
    }
  }

  //first GravMass
  theIO->GravMassBuffer[0][0] = gravMass_r(theGravMasses,0);
  theIO->GravMassBuffer[0][1] = gravMass_phi(theGravMasses,0);
  theIO->GravMassBuffer[0][2] = gravMass_M(theGravMasses,0);
  theIO->GravMassBuffer[0][3] = gravMass_omega(theGravMasses,0);
  theIO->GravMassBuffer[0][4] = gravMass_E(theGravMasses,0);
  theIO->GravMassBuffer[0][5] = gravMass_L(theGravMasses,0);
  theIO->GravMassBuffer[0][6] = gravMass_vr(theGravMasses,0);
  theIO->GravMassBuffer[0][7] = gravMass_Fr(theGravMasses,0);
  theIO->GravMassBuffer[0][8] = gravMass_Fp(theGravMasses,0);
  //second GravMass
  theIO->GravMassBuffer[1][0] = gravMass_r(theGravMasses,1);
  theIO->GravMassBuffer[1][1] = gravMass_phi(theGravMasses,1);
  theIO->GravMassBuffer[1][2] = gravMass_M(theGravMasses,1);
  theIO->GravMassBuffer[1][3] = gravMass_omega(theGravMasses,1);
  theIO->GravMassBuffer[1][4] = gravMass_E(theGravMasses,1);
  theIO->GravMassBuffer[1][5] = gravMass_L(theGravMasses,1);
  theIO->GravMassBuffer[1][6] = gravMass_vr(theGravMasses,1);
  theIO->GravMassBuffer[1][7] = gravMass_Fr(theGravMasses,1);
  theIO->GravMassBuffer[1][8] = gravMass_Fp(theGravMasses,1);

}

void io_readbuf(struct IO *theIO,struct Cell ***theCells,struct Sim *theSim,struct GravMass * theGravMasses){
  int NUM_Q = sim_NUM_Q(theSim);
  int i,j,k,q;
  int index=0;
  for (k=0; k<sim_N(theSim,Z_DIR); k++) {
    for (i=0; i<sim_N(theSim,R_DIR);i++){
      for (j=0; j<sim_N_p(theSim,i);j++){
        cell_set_tiph(theCells,i,j,k,theIO->buffer[index][0]);
        for (q=0;q<NUM_Q;q++){
          cell_set_prim(theCells,i,j,k,q,theIO->buffer[index][q+3]);
        }
        index += 1;
      }
    }
  }
  //first GravMass
  int p = 0;
  gravMass_set_chkpt(theGravMasses,p,theIO->GravMassBuffer[p][0],theIO->GravMassBuffer[p][1],theIO->GravMassBuffer[p][2],theIO->GravMassBuffer[p][3], theIO->GravMassBuffer[p][4],theIO->GravMassBuffer[p][5],theIO->GravMassBuffer[p][6],theIO->GravMassBuffer[p][7],theIO->GravMassBuffer[p][8]);
  //gravMass_set_chkpt(theGravMasses,p,theIO->GravMassBuffer[p][0],theIO->GravMassBuffer[p][1],theIO->GravMassBuffer[p][2],theIO->GravMassBuffer[p][3]);
	//second GravMass
  p=1;
  gravMass_set_chkpt(theGravMasses,p,theIO->GravMassBuffer[p][0],theIO->GravMassBuffer[p][1],theIO->GravMassBuffer[p][2],theIO->GravMassBuffer[p][3], theIO->GravMassBuffer[p][4],theIO->GravMassBuffer[p][5],theIO->GravMassBuffer[p][6],theIO->GravMassBuffer[p][7],theIO->GravMassBuffer[p][8]);
  //gravMass_set_chkpt(theGravMasses,p,theIO->GravMassBuffer[p][0],theIO->GravMassBuffer[p][1],theIO->GravMassBuffer[p][2],theIO->GravMassBuffer[p][3]);

}
