#define IO_PRIVATE_DEFS
#define DATASETNAME 	"DoubleArray" 
#define RANK   2

#include <string.h>
#include <stdlib.h>
#include "../Headers/Sim.h"
#include "../Headers/Cell.h"
#include "hdf5.h"
#include "../Headers/IO.h"
#include "../Headers/header.h"

struct IO *io_create(struct Sim *theSim) {
  struct IO *theIO = (struct IO *) malloc(sizeof(struct IO));

  return(theIO);
}

void io_destroy(struct IO *theIO){
  if(theIO->buffer != NULL) {
    free(theIO->buffer[0]);
    free(theIO->buffer);
  }
  free(theIO);
}
