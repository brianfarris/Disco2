#define PLANET_PRIVATE_DEFS
#include <stdlib.h>
#include <math.h>
#include "../Headers/GravMass.h"
#include "../Headers/header.h"


void GravMass_set_Fr( struct GravMass * theGravMasses, int p, double FFr){
  theGravMasses[p].Fr = FFr;
}

void GravMass_set_Fp( struct GravMass * theGravMasses, int p, double FFp){
  theGravMasses[p].Fp = FFp;
}

void GravMass_set_omega( struct GravMass * theGravMasses, int p, double Om){
  theGravMasses[p].omega = Om;
}
void GravMass_set_Ltot( struct GravMass * theGravMasses, int p, double Ltot){
  theGravMasses[p].Ltot = Ltot;
}
