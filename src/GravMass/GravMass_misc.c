#define PLANET_PRIVATE_DEFS
#include <stdlib.h>
#include <math.h>
#include "../Headers/GravMass.h"
#include "../Headers/header.h"

void gravMass_clean_pi(struct GravMass * theGravMasses){
   int p;
   for( p=0 ; p<NP ; ++p ){
      double phi = theGravMasses[p].phi;
      while( phi > 2.*M_PI ) phi -= 2.*M_PI;
      while( phi < 0.0 ) phi += 2.*M_PI;
      theGravMasses[p].phi = phi;
   }
}

void gravMass_copy(struct GravMass * theGravMasses){
  int p;
  for( p=0 ; p<NP ; ++p ){
    theGravMasses[p].RK_r   = theGravMasses[p].r;
    theGravMasses[p].RK_phi = theGravMasses[p].phi;
    theGravMasses[p].RK_M   = theGravMasses[p].M;
    theGravMasses[p].RK_E   = theGravMasses[p].E;
    theGravMasses[p].RK_L   = theGravMasses[p].L;
    theGravMasses[p].RK_vr  = theGravMasses[p].vr;
  }
}


