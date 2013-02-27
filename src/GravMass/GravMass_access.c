#define PLANET_PRIVATE_DEFS
#include <stdlib.h>
#include <math.h>
#include "../Headers/GravMass.h"
#include "../Headers/header.h"

double gravMass_r(struct GravMass * theGravMasses,int p){
  return(theGravMasses[p].r);
}
double gravMass_phi(struct GravMass * theGravMasses,int p){
  return(theGravMasses[p].phi);
}
double gravMass_M(struct GravMass * theGravMasses,int p){
  return(theGravMasses[p].M);
}