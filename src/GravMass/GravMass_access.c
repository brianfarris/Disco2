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
double gravMass_dist(struct GravMass * theGravMasses,int p,double r, double phi, double z){
    double xbh = theGravMasses[p].r * cos(theGravMasses[p].phi);
    double ybh = theGravMasses[p].r * sin(theGravMasses[p].phi);
    double x = r * cos(phi);
    double y = r * sin(phi);
    return(sqrt((x-xbh)*(x-xbh) + (y-ybh)*(y-ybh) + z*z));
}
