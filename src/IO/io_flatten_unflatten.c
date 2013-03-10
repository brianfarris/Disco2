#define IO_PRIVATE_DEFS
#define H5FILE_NAME     "testout.h5"
#define DATASETNAME 	"DoubleArray" 
#define RANK   2

#include <string.h>
#include <stdlib.h>
#include "../Headers/Sim.h"
#include "../Headers/Cell.h"
#include "hdf5.h"
#include "../Headers/IO.h"
#include "../Headers/header.h"

void io_flattened_prim(struct IO *io_pointer,struct Cell ***theCells,struct Sim *theSim){
  int NUM_Q = sim_NUM_Q(theSim);
  int i,j,k,q;
  int index=0;
  for (k=sim_Nghost_min(theSim,Z_DIR); k<(sim_N(theSim,Z_DIR)-sim_Nghost_max(theSim,Z_DIR)); k++) {
    for (i=sim_Nghost_min(theSim,R_DIR); i<(sim_N(theSim,R_DIR)-sim_Nghost_max(theSim,R_DIR));i++){
      for (j=0; j<sim_N_p(theSim,i);j++){
           io_pointer->primitives[index][0] = cell_tiph(cell_single(theCells,i,j,k));
          io_pointer->primitives[index][1] = 0.5*(sim_FacePos(theSim,i-1,R_DIR)+sim_FacePos(theSim,i,R_DIR));
          io_pointer->primitives[index][2] = 0.5*(sim_FacePos(theSim,k-1,Z_DIR)+sim_FacePos(theSim,k,Z_DIR));
        for (q=0;q<NUM_Q;q++){
          io_pointer->primitives[index][q+3] = cell_prim(cell_single(theCells,i,j,k),q); 
        }
        index += 1;
      }
    }
  }
}

void io_unflattened_prim(struct IO *io_pointer,struct Cell ***theCells,struct Sim *theSim){
  int NUM_Q = sim_NUM_Q(theSim);
  int i,j,k,q;
  int index=0;
  for (k=sim_Nghost_min(theSim,Z_DIR); k<(sim_N(theSim,Z_DIR)-sim_Nghost_max(theSim,Z_DIR)); k++) {
    for (i=sim_Nghost_min(theSim,R_DIR); i<(sim_N(theSim,R_DIR)-sim_Nghost_max(theSim,R_DIR));i++){
      for (j=0; j<sim_N_p(theSim,i);j++){
        cell_set_tiph(theCells,i,j,k,io_pointer->primitives[index][0]);
        for (q=0;q<NUM_Q;q++){
          cell_set_prim(theCells,i,j,k,q,io_pointer->primitives[index][q+3]);
        }
        index += 1;
      }
    }
  }
}
