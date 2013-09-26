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
  } else if  (sim_GravMassType(theSim)==LIVEBINARY){
	return(&gravMass_init_livebinary);
  } else{
    printf("ERROR - this is not the gravMass you are looking for! (look at GravMass_init.c)\n");
    exit(0);
  }
}

void gravMass_init_none(struct GravMass * theGravMasses,struct Sim * theSim){
  theGravMasses[0].M   = 0.0;
  theGravMasses[0].r   = 0.0;
  theGravMasses[0].phi = 0.0;
  theGravMasses[0].omega = 0.0;
  //theGravMasses[0].E = 0.0;
  //theGravMasses[0].L = 0.0;
  //theGravMasses[0].vr = 0.0; 
  //theGravMasses[0].Fr = 0.0; // do I need this here DD?
  ///theGravMasses[0].Fp = 0.0;  	
	
  theGravMasses[1].M   = 0.0;
  theGravMasses[1].r   = 0.0;
  theGravMasses[1].phi = 0.0;
  theGravMasses[1].omega = 0.0;
  //theGravMasses[1].E = 0.0;
  //theGravMasses[1].L = 0.0;
  //theGravMasses[1].vr = 0.0;
  //theGravMasses[1].Fr = 0.0; // do I need this here DD?
 // theGravMasses[1].Fp = 0.0;

}

void gravMass_init_single(struct GravMass * theGravMasses, struct Sim * theSim){
  theGravMasses[0].M   = 1.0;
  theGravMasses[0].r   = 0.0;
  theGravMasses[0].phi = 0.0;
  theGravMasses[0].omega = 0.0;
  //theGravMasses[0].E = 0.0;
  //theGravMasses[0].L = 0.0;
  //theGravMasses[0].vr = 0.0; 
  //theGravMasses[0].Fr = 0.0; // do I need this here DD?
  //theGravMasses[0].Fp = 0.0;  	
	
  theGravMasses[1].M   = 0.0;
  theGravMasses[1].r   = 0.0;
  theGravMasses[1].phi = 0.0;
  theGravMasses[1].omega = 0.0;
  //theGravMasses[1].E = 0.0;
  //theGravMasses[1].L = 0.0; 
  //theGravMasses[1].vr = 0.0;
  //theGravMasses[1].Fr = 0.0; // do I need this here DD?
  //theGravMasses[1].Fp = 0.0;
	

}

void gravMass_init_binary(struct GravMass * theGravMasses,struct Sim * theSim){
  double Mtotal = 1.0;
  double sep = sim_sep0(theSim);
  //double massratio = 1.0;
  double massratio = sim_MassRatio(theSim);
  double M0 = Mtotal/(1.+massratio);
  double M1 = Mtotal/(1.+1./massratio);
  double r0 = M1/Mtotal*sep;
  double r1 = M0/Mtotal*sep;
  double om2 = 1./sep/sep/sep;
	
  double L0  = M0*r0*r0*sqrt(om2);
  double L1 = M1*r1*r1*sqrt(om2);
	
  double v0 = 0.0;
  double v1 = v0*M1/M0;
  double E  = .5*M0*v0*v0 + .5*M1*v1*v1 + .5*L0*L0/M0/r0/r0 + .5*L1*L1/M1/r1/r1 - M0*M1/(r0+r1);

  theGravMasses[0].M   = M0;
  theGravMasses[0].r   = r0;
  theGravMasses[0].phi = 0.0;
  theGravMasses[0].omega = sqrt(om2);
  theGravMasses[0].E = E;
  theGravMasses[0].L = L0;
  theGravMasses[0].vr = 0.0; 
  //theGravMasses[0].Fr = 0.0; // do I need this here DD?
  //theGravMasses[0].Fp = 0.0;
	
  theGravMasses[1].M   = M1;
  theGravMasses[1].r   = r1;
  theGravMasses[1].phi = M_PI;
  theGravMasses[1].omega = sqrt(om2);
  theGravMasses[1].E = E;
  theGravMasses[1].L = L1;
  theGravMasses[1].vr = 0.0;
  //theGravMasses[1].Fr = 0.0; // do I need this here DD?
  //theGravMasses[1].Fp = 0.0;
}

// Add plant radial velocity, energy and angular momentum to GravMass structure so binary paramaters can be updated
void gravMass_init_livebinary(struct GravMass * theGravMasses,struct Sim * theSim){
	double Mtotal = 1.0;
	double sep = sim_sep0(theSim);
	double massratio = sim_MassRatio(theSim);
	double M0 = Mtotal/(1.+massratio);
	double M1 = Mtotal/(1.+1./massratio);
	double r0 = M1/Mtotal*sep;
	double r1 = M0/Mtotal*sep;
	double om = sqrt(1./sep/sep/sep);
	
	double L0 = M0*r0*r0*om;
	double L1 = M1*r1*r1*om;
	
	double v0 = 0.0;
	double v1 = v0*M1/M0;
	//double E  = .5*M0*v0*v0 + .5*M1*v1*v1 + .5*L0*L0/M0/r0/r0 + .5*L1*L1/M1/r1/r1 - M0*M1/(r0+r1); //baryocentric 
	double E  = -0.5*M0*M1/(r0+r1); //barycentric 
	//double E  = .5*M0*v0*v0 + .5*M1*v1*v1 + .5*L0*L0/M0/r0/r0 + .5*L1*L1/M1/r1/r1 - (M0+M1)/(r0+r1); // see murray dermott
	
	theGravMasses[0].M   = M0;
	theGravMasses[0].r   = r0;
	theGravMasses[0].phi = 0.0;
	theGravMasses[0].omega = om;
	theGravMasses[0].E = E;
	theGravMasses[0].L = L0;
	theGravMasses[0].vr = v0;
	theGravMasses[0].Fr = 0.0; 
	theGravMasses[0].Fp = 0.0;
	
	theGravMasses[1].M   = M1;
	theGravMasses[1].r   = r1;
	theGravMasses[1].phi = M_PI;
	theGravMasses[1].omega = om;
	theGravMasses[1].E = E;
	theGravMasses[1].L = L1;
	theGravMasses[1].vr = v1;
	theGravMasses[1].Fr = 0.0; 
	theGravMasses[1].Fp = 0.0;
}



 void gravMass_set_chkpt(struct GravMass * theGravMasses,int p,double r,double phi,double M, double omega, double E, double L, double vr, double Fr, double Fp){
	 theGravMasses[p].r = r; 
	 theGravMasses[p].phi = phi;
	 theGravMasses[p].M = M;
	 theGravMasses[p].omega = omega;
	 theGravMasses[p].E = E;
	 theGravMasses[p].L = L;
	 theGravMasses[p].vr = vr; 
	 theGravMasses[p].Fr = Fr; 
	 theGravMasses[p].Fp = Fp;  		
 }



/*
void gravMass_set_chkpt(struct GravMass * theGravMasses,int p,double r,double phi,double M, double omega){
	theGravMasses[p].r = r; 
	theGravMasses[p].phi = phi;
	theGravMasses[p].M = M;
	theGravMasses[p].omega = omega; 		
}
 */
