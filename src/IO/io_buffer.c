#define IO_PRIVATE_DEFS
#include <string.h>
#include <stdlib.h>
#include "../Headers/Sim.h"
#include "../Headers/Cell.h"
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

void io_setbuf(struct IO *theIO,struct Cell ***theCells,struct Sim *theSim){
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
}

void io_readbuf(struct IO *theIO,struct Cell ***theCells,struct Sim *theSim){
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
}
