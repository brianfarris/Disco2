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
double gravMass_Omega_eff(double r, double phi,int NumGravMass, struct GravMass * theGravMasses){
    double Omega2_eff = 0.0;
    int p;
    for (p=0;p<NumGravMass; ++p ){
        double m = gravMass_M(theGravMasses,p);
        double dist_bh = gravMass_dist(theGravMasses,p,r,phi,0.);
        Omega2_eff += m/pow(dist_bh,3);
    }
    return(sqrt(Omega2_eff));
}
