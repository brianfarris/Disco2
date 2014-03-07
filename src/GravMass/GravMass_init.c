#define PLANET_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/GravMass.h"
#include "../Headers/Sim.h"
#include "../Headers/header.h"

void (*gravMass_init_ptr(struct Sim * theSim))(struct GravMass *,struct Sim *){
  if (sim_GravMassType(theSim)==NONE){
    return(&gravMass_init_none);
  } else if (sim_GravMassType(theSim)==SINGLE){
    return(&gravMass_init_single);
  } else if  (sim_GravMassType(theSim)==BINARY){
    return(&gravMass_init_binary);
  } else{
    printf("ERROR\n");
    exit(0);
  }
}

void gravMass_init_none(struct GravMass * theGravMasses,struct Sim * theSim){
  theGravMasses[0].M   = 0.0;
  theGravMasses[0].r   = 0.0;
  theGravMasses[0].phi = 0.0;
  theGravMasses[0].omega = 0.0;
  theGravMasses[0].Mdot = 0.0;
  theGravMasses[1].M   = 0.0;
  theGravMasses[1].r   = 0.0;
  theGravMasses[1].phi = 0.0;
  theGravMasses[1].omega = 0.0;
  theGravMasses[1].Mdot = 0.0;
}

void gravMass_init_single(struct GravMass * theGravMasses,struct Sim * theSim){
  theGravMasses[0].M   = 1.0;
  theGravMasses[0].r   = 0.0;
  theGravMasses[0].phi = 0.0;
  theGravMasses[0].omega = 0.0;
  theGravMasses[0].Mdot = 0.0;
  theGravMasses[1].M   = 0.0;
  theGravMasses[1].r   = 1.0;
  theGravMasses[1].phi = 0.0;
  theGravMasses[1].omega = 0.0;
  theGravMasses[1].Mdot = 0.0;

}

void gravMass_init_binary(struct GravMass * theGravMasses,struct Sim * theSim){
  double Mtotal = 1.0;
  double sep = 1.0;
  //double massratio = 1.0;
  double massratio = sim_MassRatio(theSim);
  theGravMasses[0].OrbShrinkTscale = sim_OrbShrinkTscale(theSim);
  theGravMasses[1].OrbShrinkTscale = sim_OrbShrinkTscale(theSim);
  theGravMasses[0].OrbShrinkT0 = sim_OrbShrinkT0(theSim);
  theGravMasses[1].OrbShrinkT0 = sim_OrbShrinkT0(theSim);
  double M0 = Mtotal/(1.+massratio);
  double M1 = Mtotal/(1.+1./massratio);
  double r0 = M1/Mtotal*sep;
  double r1 = M0/Mtotal*sep;
  double om2 = 1./sep/sep/sep;

  theGravMasses[0].M   = M0;
  theGravMasses[0].r   = r0;
  theGravMasses[0].phi = 0.0;
  theGravMasses[0].omega = sqrt(om2);
  theGravMasses[0].Mdot = 0.0;
  theGravMasses[0].OrbShrinkTscale = sim_OrbShrinkTscale(theSim);
  theGravMasses[0].OrbShrinkT0 = sim_OrbShrinkT0(theSim);
  theGravMasses[1].M   = M1;
  theGravMasses[1].r   = r1;
  theGravMasses[1].phi = M_PI;
  theGravMasses[1].omega = sqrt(om2);
  theGravMasses[1].Mdot = 0.0;
  theGravMasses[1].OrbShrinkTscale = sim_OrbShrinkTscale(theSim);
  theGravMasses[1].OrbShrinkT0 = sim_OrbShrinkT0(theSim);
}

void gravMass_set_chkpt(struct GravMass * theGravMasses,int p,double r,double phi,double M, double omega){
  theGravMasses[p].r = r; 
  theGravMasses[p].phi = phi;
  theGravMasses[p].M = M;
  theGravMasses[p].omega = omega;
}
