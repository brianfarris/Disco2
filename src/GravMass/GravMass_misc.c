#define PLANET_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/GravMass.h"
#include "../Headers/Sim.h"
#include "../Headers/header.h"

void gravMass_clean_pi(struct GravMass * theGravMasses,struct Sim * theSim){
   int p;
   for( p=0 ; p<sim_NumGravMass(theSim) ; ++p ){
      double phi = theGravMasses[p].phi;
      while( phi > 2.*M_PI ) phi -= 2.*M_PI;
      while( phi < 0.0 ) phi += 2.*M_PI;
      theGravMasses[p].phi = phi;
   }
}


void gravMass_copy(struct GravMass * theGravMasses,struct Sim * theSim){
  int p;
  for( p=0 ; p<sim_NumGravMass(theSim) ; ++p ){
    theGravMasses[p].RK_r   = theGravMasses[p].r;
    theGravMasses[p].RK_phi = theGravMasses[p].phi;
    theGravMasses[p].RK_M   = theGravMasses[p].M;
    theGravMasses[p].RK_omega   = theGravMasses[p].omega;
//    theGravMasses[p].RK_E   = theGravMasses[p].E;
//    theGravMasses[p].RK_L   = theGravMasses[p].L;
//    theGravMasses[p].RK_vr  = theGravMasses[p].vr;
  }
}


void gravMass_move(struct GravMass * theGravMasses,double t,double dt){
  double m0 = theGravMasses[0].M;
  double m1 = theGravMasses[1].M;
  double M = m0+m1;
  double mu = m0*m1/M; 
  double a_0 = 1.0;
  double tau_0 = 1000.0;//5./256. * pow(a_0,4)/(mu * M*M);
  double a = a_0 * pow((1.0 - (t-2300.)/tau_0),0.25);
  double omega = pow(a/M,-1.5);
  //printf("t-2300: %e, a: %e, omega: %e\n",t-2300.,a,omega);
  /*
     theGravMasses[0].phi += theGravMasses[0].omega*dt;
     theGravMasses[1].phi += theGravMasses[1].omega*dt;
     */
  theGravMasses[0].phi += omega*dt;
  theGravMasses[1].phi += omega*dt;
}

void gravMass_update_RK( struct GravMass * theGravMasses,struct Sim * theSim, double RK){
  int i;
  for( i=0 ; i<sim_NumGravMass(theSim) ; ++i ){
    while( theGravMasses[i].phi-theGravMasses[i].RK_phi >  M_PI ) theGravMasses[i].RK_phi += 2.*M_PI;
    while( theGravMasses[i].phi-theGravMasses[i].RK_phi < -M_PI ) theGravMasses[i].RK_phi -= 2.*M_PI;

    theGravMasses[i].r   = (1.0-RK)*theGravMasses[i].r   + RK*theGravMasses[i].RK_r;
    theGravMasses[i].phi = (1.0-RK)*theGravMasses[i].phi + RK*theGravMasses[i].RK_phi;
    theGravMasses[i].M   = (1.0-RK)*theGravMasses[i].M   + RK*theGravMasses[i].RK_M;
    //theGravMasses[i].L   = (1.0-RK)*theGravMasses[i].L   + RK*theGravMasses[i].RK_L;
    //theGravMasses[i].E   = (1.0-RK)*theGravMasses[i].E   + RK*theGravMasses[i].RK_E;
    //theGravMasses[i].vr  = (1.0-RK)*theGravMasses[i].vr  + RK*theGravMasses[i].RK_vr;
    theGravMasses[i].omega = (1.0-RK)*theGravMasses[i].omega  + RK*theGravMasses[i].RK_omega;
  }

}

void gravMass_set_Mdot( struct GravMass * theGravMasses, double Mdot, int p){
  theGravMasses[p].Mdot = Mdot;
}

void gravMass_set_total_torque( struct GravMass * theGravMasses, double total_torque, int p){
  theGravMasses[p].total_torque = total_torque;
}
