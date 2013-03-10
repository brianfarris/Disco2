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
  for (k=sim_Nghost_zmin(theSim); k<(sim_N_z(theSim)-sim_Nghost_zmax(theSim)); k++) {
    for (i=sim_Nghost_rmin(theSim); i<(sim_N_r(theSim)-sim_Nghost_rmax(theSim));i++){
      for (j=0; j<sim_N_p(theSim,i);j++){
           io_pointer->primitives[index][0] = cell_tiph(cell_single(theCells,i,j,k));
          io_pointer->primitives[index][1] = 0.5*(sim_r_faces(theSim,i-1)+sim_r_faces(theSim,i));
          io_pointer->primitives[index][2] = 0.5*(sim_z_faces(theSim,k-1)+sim_z_faces(theSim,k));
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
  for (k=sim_Nghost_zmin(theSim); k<(sim_N_z(theSim)-sim_Nghost_zmax(theSim)); k++) {
    for (i=sim_Nghost_rmin(theSim); i<(sim_N_r(theSim)-sim_Nghost_rmax(theSim));i++){
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
