#define PLANET_PRIVATE_DEFS
#include <stdlib.h>
#include <math.h>
#include "../Headers/GravMass.h"
#include "../Headers/header.h"

struct GravMass *gravMass_create(int num_gravMasses){
  struct GravMass *theGravMasses = (struct GravMass *) malloc(sizeof(struct GravMass)*num_gravMasses);
  return(theGravMasses);
}

void gravMass_destroy(struct GravMass * theGravMasses){
  free(theGravMasses);
}

