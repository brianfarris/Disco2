#define IO_PRIVATE_DEFS
#define DATASETNAME 	"DoubleArray" 
#define RANK   2

#include <string.h>
#include <stdlib.h>
#include "../Headers/Grid.h"
#include "../Headers/Cell.h"
#include "hdf5.h"
#include "../Headers/IO.h"
#include "../Headers/header.h"

struct IO *io_create(struct Grid *theGrid) {
  struct IO *theIO = (struct IO *) malloc(sizeof(struct IO));
  int NUM_Q = grid_NUM_Q(theGrid);

  int Ncells = grid_Ncells(theGrid)*(NUM_Q+3);
  int N_r = grid_N_r(theGrid);
  int N_z = grid_N_z(theGrid);

  theIO->primitives = malloc(sizeof(double *)*grid_Ncells(theGrid));
  theIO->primitives[0] = malloc(sizeof(double )*grid_Ncells(theGrid)*(NUM_Q+3));
  int i;
  for (i=0;i<grid_Ncells(theGrid);++i){
    theIO->primitives[i]= theIO->primitives[0]+i*(NUM_Q+3);
  }
  return(theIO);
}

void io_destroy(struct IO *theIO,struct Grid *theGrid){
  int i;
  free(theIO->primitives[0]);
  free(theIO->primitives);
  free(theIO);
}
