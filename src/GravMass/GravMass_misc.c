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
    theGravMasses[p].RK_E   = theGravMasses[p].E;
    theGravMasses[p].RK_L   = theGravMasses[p].L;
    theGravMasses[p].RK_Ltot   = theGravMasses[p].Ltot;
    theGravMasses[p].RK_vr  = theGravMasses[p].vr;
  }
}



// Analytic uodate via Keplers Laws - now jsut e=0 circular orbits
void gravMass_adv_anly( struct GravMass * thisGravMass , double Mp , double dt ){
  double Ms = thisGravMass->M;  //The mass of the secondary
  //-double mt = Mp/pow(1.+Ms/Mp,2.);            // ??  
  //-double E = thisGravMass->E;   // Etot                             
  double Ltot = thisGravMass->Ltot;           // Lsec/Msec
  //-double Lp = Ls*Ms/Mp; // Omega the same right? 3 1/q's to convert 1M and 2 r's works out to one Ms/Mp
  //-double Ltot = Ls; //FOR NOW Ls is tot L - DD April 7 2014 *(1.+Ms/Mp);
  //-double r1 = thisGravMass->r;
  //-double phi = thisGravMass->phi;
  //-double vr1 = thisGravMass->vr;

  //- double rt = r1 + r1*Ms/Mp;

  //either update with Delta L
  double a = Ltot*Ltot/(Ms*Ms*Mp*Mp)*(Ms+Mp);

  // Or update with Delta E                                
  //double a = -.5*Ms*Mp/Enew;               // E = -G(M0+M1)/(2a) for elliptical orbit
  // for e=0 only
  //-double e=0.0;
  //double cosE0 = 0.0;              
  //-double cosE0 = (1.-r1/a)/e;    // from r=a(1-ecos(E))        
  //-if( e <= 0.0 ) cosE0 = 0.0;   // circular orbit or imaginary e -> something weird? 
  //-double sinE0 = sqrt(fabs(1.-cosE0*cosE0));
  //-if( vr1 < 0.0 ) sinE0 = -sinE0; //???



  //-double E0 = atan2( sinE0 , cosE0 );
  //-if( E>0.0 ) E0 = log( sinE0 + cosE0 );

  //-double phi0 = 2.*atan2( sqrt(1.+e)*sin(E0/2.) , sqrt(1.-e)*cos(E0/2.) );   //standard Newtonian formulae for phi
  //if( E>0.0 ) phi0 = 2.*atan2( sqrt(1.+e)*sinh(E0/2.) , sqrt(fabs(1.-e))*cosh(E0/2.) ); //Paul
  //-if( E>0.0 ) phi0 = 2.*atan2( sqrt(e+1.)*sinh(E0/2.) , sqrt(e-1.)*cosh(E0/2.) );//above unbound; why use fabs? - if unbound then e>=1 -Dan                        
 
  //-double M0 = E0 - e*sin(E0);   // Newtonian formulae for t(=M here) in terms of Eccentric anomaly E and eccentricity e
  //-if( E>0.0 ) M0 = E0 - e*sinh(E0); // above unbound          

  //-double M = M0 + dt*sqrt(mt)/pow(a,1.5);  // advance t
  //-if( E>0.0 ) M = M0 - dt*sqrt(mt)/pow(-a,1.5);

  //-double E1 = M;
  //-double df = M-E1+e*sin(E1);
  //-if( E>0.0 ) df = M-E1+e*sinh(E1);
  //-while( fabs(df) > 1e-12 ){
  //-  double dfdE = -1.+e*cos(E1);
  //-  if( E>0.0 ) dfdE = -1.+e*cosh(E1);
  //-  double dE = -df/dfdE;
  //-  E1 += dE;
  //-  df = M-E1+e*sin(E1);
  //-  if( E>0.0 ) df = M-E1+e*sinh(E1);
  //-}

  //-rt = a*(1.-e*cos(E1));
  //-if( E>0.0 ) rt = a*(1.-e*cosh(E1));
  //-double phi1 = 2.*atan2( sqrt(1.+e)*sin(E1/2.) , sqrt(1.-e)*cos(E1/2.) );
  //-if( E>0.0 )  phi1 = 2.*atan2( sqrt(1.+e)*sinh(E1/2.) , sqrt(fabs(1.-e))*cosh(E1/2.) );


  //-thisGravMass->r = rt/(1.+Ms/Mp);
  //-double vr = sqrt((Ms+Mp)*(2./rt - 1./a) - (Ms+Mp)/a/a/a * rt*rt); //0 for circ orbit

  //-if((Ms+Mp)*(2./rt - 1./a) - (Ms+Mp)/a/a/a * rt*rt < 0.0 ) vr = 0.0;
  //-if( sin(E1) < 0.0 ) vr = -vr; //???
  //*vrt = vr2;
  //-thisGravMass->vr = vr/(1.+Ms/Mp);
  thisGravMass->r = a/(1.+Ms/Mp);
  thisGravMass->omega = pow(a,(-3./2.));           //Is this a bad idea?
}





//Analytic Kepler    
void gravMass_move( struct Sim * theSim, struct GravMass * theGravMasses, double t, double dt ){
  if (sim_GravMassType(theSim)==LIVEBINARY){
    //printf("LIVEBINARY");
    double Mp = theGravMasses[0].M; //primary
    double Ms = theGravMasses[1].M; //secondary

    //Here we need to analytically update a, phi, vr from Kepler's equations 
    gravMass_adv_anly(&(theGravMasses[1]), Mp, dt);
    // Or just push the binary together arbitrarily for testing
    //gravMass_adv_arb(&(theGravMasses[1]), Mp, dt); //TESTING                                                                        
    theGravMasses[0].r   = theGravMasses[1].r*Ms/Mp;
    theGravMasses[0].vr  = theGravMasses[1].vr*Ms/Mp;
    theGravMasses[0].omega  = theGravMasses[1].omega;

    theGravMasses[1].phi += theGravMasses[1].omega*dt; //Instead of updating in adv_anly ONLY GOOD FOR CIRC?
    theGravMasses[0].phi = theGravMasses[1].phi + M_PI;
 }else{
    theGravMasses[1].phi += theGravMasses[1].omega*dt;
    //theGravMasses[1].phi += theGravMasses[1].omega*dt               
    theGravMasses[0].phi = theGravMasses[1].phi + M_PI;
  }
}




/*
void gravMass_move(struct Sim * theSim,struct GravMass * theGravMasses,double t,double dt){
  double m0 = theGravMasses[0].M;
  double m1 = theGravMasses[1].M;
  double M = m0+m1;
  double mu = m0*m1/M; 
  double a_0 = 1.0;
  double OrbShrinkTscale = sim_OrbShrinkTscale(theSim);
  double OrbShrinkT0 = sim_OrbShrinkT0(theSim);
  double tau_0 = OrbShrinkTscale;
  double a; 
  double omega;
  if ((t-OrbShrinkT0)/tau_0<1.){
    a = a_0 * pow((1.0 - (t-OrbShrinkT0)/tau_0),0.25);
    omega = pow(a/M,-1.5);
  } else{
    a = 0.0;
    omega = 0.0;
  }

  //printf("a: %e, M: %e, omega: %e\n",a,M,omega);
  double r0,r1;
  if (M<1.e-12){
    r0 = 0.0;
    r1=0.0;
  } else{
    r0 = m1/M*a;
    r1 = m0/M*a;
  }
     /*
     theGravMasses[0].phi += theGravMasses[0].omega*dt;
     theGravMasses[1].phi += theGravMasses[1].omega*dt;
     */ /*
  theGravMasses[0].phi += omega*dt;
  theGravMasses[1].phi += omega*dt;
  theGravMasses[0].r = r0;
  theGravMasses[1].r = r1;
  theGravMasses[0].omega = omega;
  theGravMasses[1].omega = omega;
}
*/


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
    ///Set a density scale for feedback to holes                                 
    //    double Rcut = sim_Rcut(theSim);
    // double dens_scale = ( sim_Mdsk_o_Ms(theSim)/( 1.+1./sim_MassRatio(theSim) ) )/( M_PI*sim_sep0(theSim)*sim_sep0(theSim) );
    //if (   time_global <=    2.*M_PI*(sim_tmig_on(theSim) + sim_tramp(theSim))   ){
    //  dens_scale *= (time_global - 2*M_PI*sim_tmig_on(theSim))/( sim_tramp(theSim)*2*M_PI );
    //  // ABOVE^ Allow migration to turn on slowly 
    //}
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
