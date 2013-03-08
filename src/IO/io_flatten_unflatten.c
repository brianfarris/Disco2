#define IO_PRIVATE_DEFS
#define H5FILE_NAME     "testout.h5"
#define DATASETNAME 	"DoubleArray" 
#define RANK   2

#include <string.h>
#include <stdlib.h>
#include "../Headers/Grid.h"
#include "../Headers/Cell.h"
#include "hdf5.h"
#include "../Headers/IO.h"
#include "../Headers/header.h"

void io_flattened_prim(struct IO *io_pointer,struct Cell ***theCells,struct Grid *theGrid){
  int NUM_Q = grid_NUM_Q(theGrid);
  int i,j,k,q;
  int index=0;
  for (k=grid_Nghost_zmin(theGrid); k<(grid_N_z(theGrid)-grid_Nghost_zmax(theGrid)); k++) {
    for (i=grid_Nghost_rmin(theGrid); i<(grid_N_r(theGrid)-grid_Nghost_rmax(theGrid));i++){
      for (j=0; j<grid_N_p(theGrid,i);j++){
           io_pointer->primitives[index][0] = cell_tiph(cell_single(theCells,i,j,k));
          io_pointer->primitives[index][1] = 0.5*(grid_r_faces(theGrid,i-1)+grid_r_faces(theGrid,i));
          io_pointer->primitives[index][2] = 0.5*(grid_z_faces(theGrid,k-1)+grid_z_faces(theGrid,k));
        for (q=0;q<NUM_Q;q++){
          io_pointer->primitives[index][q+3] = cell_prim(cell_single(theCells,i,j,k),q); 
        }
        index += 1;
      }
    }
  }
}

void io_unflattened_prim(struct IO *io_pointer,struct Cell ***theCells,struct Grid *theGrid){
  int NUM_Q = grid_NUM_Q(theGrid);
  int i,j,k,q;
  int index=0;
  for (k=grid_Nghost_zmin(theGrid); k<(grid_N_z(theGrid)-grid_Nghost_zmax(theGrid)); k++) {
    for (i=grid_Nghost_rmin(theGrid); i<(grid_N_r(theGrid)-grid_Nghost_rmax(theGrid));i++){
      for (j=0; j<grid_N_p(theGrid,i);j++){
        cell_set_tiph(theCells,i,j,k,io_pointer->primitives[index][0]);
        for (q=0;q<NUM_Q;q++){
          cell_set_prim(theCells,i,j,k,q,io_pointer->primitives[index][q+3]);
        }
        index += 1;
      }
    }
  }
}
