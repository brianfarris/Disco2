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
double gravMass_omega(struct GravMass * theGravMasses,int p){
  return(theGravMasses[p].omega);
}
double gravMass_Mdot(struct GravMass * theGravMasses,int p){
  return(theGravMasses[p].Mdot);
}
double gravMass_Macc(struct GravMass * theGravMasses,int p){
  return(theGravMasses[p].Macc);
}
double gravMass_total_torque(struct GravMass * theGravMasses,int p){
  return(theGravMasses[p].total_torque);
}
double gravMass_E(struct GravMass * theGravMasses,int p){
  return(theGravMasses[p].E);
}
double gravMass_L(struct GravMass * theGravMasses,int p){
  return(theGravMasses[p].L);
}
double gravMass_Ltot(struct GravMass * theGravMasses,int p){
  return(theGravMasses[p].Ltot);
}
double gravMass_vr(struct GravMass * theGravMasses,int p){
  return(theGravMasses[p].vr);
}
double gravMass_Fr(struct GravMass * theGravMasses,int p){
  return(theGravMasses[p].Fr);
}
double gravMass_Fp(struct GravMass * theGravMasses,int p){
  return(theGravMasses[p].Fp);
}
