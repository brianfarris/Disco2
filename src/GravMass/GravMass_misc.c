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
    theGravMasses[p].RK_r      = theGravMasses[p].r;
    theGravMasses[p].RK_phi    = theGravMasses[p].phi;
    theGravMasses[p].RK_M      = theGravMasses[p].M;
    theGravMasses[p].RK_omega  = theGravMasses[p].omega;
    theGravMasses[p].RK_E      = theGravMasses[p].E;
    theGravMasses[p].RK_L      = theGravMasses[p].L;
    theGravMasses[p].RK_Ltot   = theGravMasses[p].Ltot;
    theGravMasses[p].RK_vr     = theGravMasses[p].vr;
  }
}



// Analytic uodate via Keplers Laws - now jsut e=0 circular orbits
void gravMass_adv_anly( struct GravMass * thisGravMass , double Mp , double dt ){
    double Ms = thisGravMass->M; 
    double Ltot = thisGravMass->Ltot;
    double a = Ltot*Ltot/(Ms*Ms*Mp*Mp)*(Ms+Mp);

    thisGravMass->r = a/(1.+Ms/Mp);
    thisGravMass->omega = pow(a,(-3./2.));           //Is this a bad idea?
}

// Forced uodate
void gravMass_adv_arb( struct Sim * theSim, struct GravMass * thisGravMass , double Mp , double t, double dt ){
    double anew;
    double Ms = thisGravMass->M;
    double rs0 = 1./(1.+Ms/Mp);
    double nu = sim_EXPLICIT_VISCOSITY(theSim);
    //double drdt;
    double rate = -sim_vRate(theSim)*1.5 * nu/rs0;
    double dexp = 2./3.;
    //double aold = thisGravMass->r*(1.+Ms/Mp); 
    double tmax = sim_tmaxOrb(theSim) * 2.*M_PI;
    anew = pow( (1. + rate * (t - tmax)/dexp), dexp);
    //anew = pow( (1. + rate * (t - tmax))/dexp, dexp);

    /*
    if ( time_global > 2.*M_PI*sim_tmig_on(theSim) ){
      drdt = -sim_vRate(theSim)*1.5 * nu/rs0;
    }else{
      drdt = 0.0;
    }
    */
    //thisGravMass->r += drdt * dt;
    //double a = thisGravMass->r*(1.+Ms/Mp);
    //thisGravMass->omega = pow(a,(-3./2.));           //Is this a bad idea?                          
    thisGravMass->r = anew/(1.+Ms/Mp);
    thisGravMass->omega = pow(anew,-1.5);
}





//Analytic Kepler    
void gravMass_move( struct Sim * theSim, struct GravMass * theGravMasses, double t, double dt ){
  if (sim_GravMassType(theSim)==LIVEBINARY){
    //printf("LIVEBINARY");
    double Mp = theGravMasses[0].M; //primary
    double Ms = theGravMasses[1].M; //secondary

    //Here we need to analytically update a, phi, vr from Kepler's equations 
    //gravMass_adv_anly(&(theGravMasses[1]), Mp, dt);
    // Or just push the binary together arbitrarily for testing
    gravMass_adv_arb(theSim, &(theGravMasses[1]), Mp, t, dt); //TESTING 
                                                                       
    theGravMasses[0].r   = theGravMasses[1].r*Ms/Mp;
    theGravMasses[0].omega  = theGravMasses[1].omega;

    theGravMasses[1].phi += theGravMasses[1].omega*dt; //Instead of updating in adv_anly ONLY GOOD FOR CIRC?
    theGravMasses[0].phi = theGravMasses[1].phi + M_PI;
 }else{
    theGravMasses[1].phi += theGravMasses[1].omega*dt;
    theGravMasses[0].phi = theGravMasses[1].phi + M_PI;
  }
}





void gravMass_update_RK( struct GravMass * theGravMasses,struct Sim * theSim, double RK, double dt){
  int i;
  for( i=0 ; i<sim_NumGravMass(theSim) ; ++i ){
    while( theGravMasses[i].phi-theGravMasses[i].RK_phi >  M_PI ) theGravMasses[i].RK_phi += 2.*M_PI;
    while( theGravMasses[i].phi-theGravMasses[i].RK_phi < -M_PI ) theGravMasses[i].RK_phi -= 2.*M_PI;

    theGravMasses[i].r     = (1.0-RK)*theGravMasses[i].r      + RK*theGravMasses[i].RK_r;
    theGravMasses[i].phi   = (1.0-RK)*theGravMasses[i].phi    + RK*theGravMasses[i].RK_phi;
    theGravMasses[i].M     = (1.0-RK)*theGravMasses[i].M      + RK*theGravMasses[i].RK_M;
    theGravMasses[i].L     = (1.0-RK)*theGravMasses[i].L      + RK*theGravMasses[i].RK_L;
    theGravMasses[i].Ltot  = (1.0-RK)*theGravMasses[i].Ltot   + RK*theGravMasses[i].RK_Ltot;
    theGravMasses[i].E     = (1.0-RK)*theGravMasses[i].E      + RK*theGravMasses[i].RK_E;
    theGravMasses[i].vr    = (1.0-RK)*theGravMasses[i].vr     + RK*theGravMasses[i].RK_vr;
    theGravMasses[i].omega = (1.0-RK)*theGravMasses[i].omega  + RK*theGravMasses[i].RK_omega;
  }
  if ( (sim_GravMassType(theSim)==LIVEBINARY)  && (time_global > 2.*M_PI*sim_tmig_on(theSim)) ){

    //toal torque gives the torque on the binary due to the gas wrt to r=0 
    double Tot_Trq =  gravMass_total_torque(theGravMasses,0);
    double Lbin    =  gravMass_Ltot(theGravMasses,1); //theGravMasses[1].Ltot;
    //IMPORTANT 0 and 1 Ltot are the same but use 1 only! see adv_anyl() above    

    Lbin  += Tot_Trq*dt;      
    GravMass_set_Ltot(theGravMasses, 1, Lbin);
  }
}

void gravMass_set_Mdot( struct GravMass * theGravMasses, double Mdot, int p){
  theGravMasses[p].Mdot = Mdot;
}

void gravMass_set_Macc( struct GravMass * theGravMasses, double Macc, int p){
  theGravMasses[p].Macc = Macc;
}

void gravMass_set_total_torque( struct GravMass * theGravMasses, double total_torque, int p){
  theGravMasses[p].total_torque = total_torque;
}
