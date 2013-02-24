#define PLANET_PRIVATE_DEFS
#include <stdlib.h>
#include <math.h>
#include "../Headers/GravMass.h"
#include "../Headers/header.h"

void gravMass_initialize(struct GravMass * theGravMasses){
  theGravMasses[0].M   = 1.0;
  theGravMasses[0].E   = 0.0;
  theGravMasses[0].L   = 0.0;
  theGravMasses[0].r   = 0.0;
  theGravMasses[0].phi = 0.0;
}

