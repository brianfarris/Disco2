#define IO_PRIVATE_DEFS
#define H5FILE_NAME     "testout.h5"
#define DATASETNAME 	"DoubleArray" 
#define RANK   2

#include <string.h>
#include <stdlib.h>
#include "../Headers/grid.h"
#include "../Headers/cell.h"
#include "hdf5.h"
#include "../Headers/io.h"
#include "../Headers/header.h"

struct io *io_create(struct Grid *theGrid) {
  int Ncells = grid_Ncells(theGrid)*(NUM_Q+3);
  int N_r = grid_N_r(theGrid);
  int N_z = grid_N_z(theGrid);

  struct io *io_pointer = (struct io *) malloc(sizeof(struct io));

  io_pointer->primitives = malloc(sizeof(double *)*grid_Ncells(theGrid));
  io_pointer->primitives[0] = malloc(sizeof(double )*grid_Ncells(theGrid)*(NUM_Q+3));
  int i;
  for (i=0;i<grid_Ncells(theGrid);++i){
    io_pointer->primitives[i]= io_pointer->primitives[0]+i*(NUM_Q+3);
  }
  return(io_pointer);
}

void io_destroy(struct io *io_pointer,struct Grid *theGrid){
  int i;
  free(io_pointer->primitives[0]);
  free(io_pointer->primitives);
  free(io_pointer);
}
