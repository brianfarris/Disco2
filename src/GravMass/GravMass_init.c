#define PLANET_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/GravMass.h"
#include "../Headers/Grid.h"
#include "../Headers/header.h"

void (*gravMass_init_ptr(struct Grid * theGrid))(struct GravMass *){
  if (grid_GravMassType(theGrid)==NONE){
    return(&gravMass_init_none);
  } else if (grid_GravMassType(theGrid)==SINGLE){
    return(&gravMass_init_single);
  } else if  (grid_GravMassType(theGrid)==BINARY){
    return(&gravMass_init_binary);
  } else{
    printf("ERROR\n");
    exit(0);
  }
}

void gravMass_init_none(struct GravMass * theGravMasses){
  theGravMasses[0].M   = 0.0;
  theGravMasses[0].r   = 0.0;
  theGravMasses[0].phi = 0.0;
  theGravMasses[0].omega = 0.0;
  theGravMasses[1].M   = 0.0;
  theGravMasses[1].r   = 0.0;
  theGravMasses[1].phi = 0.0;
  theGravMasses[1].omega = 0.0;

}

void gravMass_init_single(struct GravMass * theGravMasses){
  theGravMasses[0].M   = 1.0;
  theGravMasses[0].r   = 0.0;
  theGravMasses[0].phi = 0.0;
  theGravMasses[0].omega = 0.0;
  theGravMasses[1].M   = 0.0;
  theGravMasses[1].r   = 0.0;
  theGravMasses[1].phi = 0.0;
  theGravMasses[1].omega = 0.0;

}

void gravMass_init_binary(struct GravMass * theGravMasses){
  double Mtotal = 1.0;
  double sep = 1.0;
  double massratio = 1.0;
  double M0 = Mtotal/(1.+massratio);
  double M1 = Mtotal/(1.+1./massratio);
  double r0 = M1/Mtotal*sep;
  double r1 = M0/Mtotal*sep;
  double om2 = 1./sep/sep/sep;

  theGravMasses[0].M   = M0;
  theGravMasses[0].r   = r0;
  theGravMasses[0].phi = 0.0;
  theGravMasses[0].omega = sqrt(om2);
  theGravMasses[1].M   = M1;
  theGravMasses[1].r   = r1;
  theGravMasses[1].phi = M_PI;
  theGravMasses[1].omega = sqrt(om2);
}

