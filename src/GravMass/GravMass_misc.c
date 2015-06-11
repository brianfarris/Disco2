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
  // Should be given the updated M, Mp, l, E. G=1
    double Ms = thisGravMass->M; 
    double Ltot = thisGravMass->Ltot;
    double Etot = thisGravMass->E;
    double asmaj = 0.5*Ms*Mp/fabs(Etot);
    //double asmaj = Ltot*Ltot/(Ms*Ms*Mp*Mp)*(Ms+Mp);
    /// Dan June 5
    double ecc = sqrt( fabs(1.- 2.*fabs(Etot)*Ltot*Ltot/pow(Ms*Mp,3)*(Ms+Mp) )  ); //E<0
    double r1 = thisGravMass->r;
    double phi = thisGravMass->phi;
    double vr1 = thisGravMass->vr;

    double rsep = r1 + r1*Ms/Mp;

    //Get inital eccentric anamoly 
    double cosE0 = (1.-rsep/asmaj)/ecc;    // from r=a(1-ecos(E))
    if( ecc <= 0.0 ) cosE0 = 0.0;   // circular orbit or imaginary e -> something weird?

    double sinE0 = sqrt(fabs(1.-cosE0*cosE0));
    if( vr1 < 0.0 ) sinE0 = -sinE0; 
    double E0 = atan2( sinE0 , cosE0 );
    if( Etot>0.0 ) E0 = log( sinE0 + cosE0 );


    //From ecc and eccentric anamoly get true anamoly
    double phi0 = 2.*atan2( sqrt(1.+ecc)*sin(E0/2.) , sqrt(1.-ecc)*cos(E0/2.) );
    //if( E>0.0 ) phi0 = 2.*atan2( sqrt(1.+e)*sinh(E0/2.) , sqrt(fabs(1.-e))*cosh(E0/2.) ); //Paul
    if( Etot>0.0 ) phi0 = 2.*atan2( sqrt(ecc+1.)*sinh(E0/2.) , sqrt(ecc-1.)*cosh(E0/2.) );//above unbound; why use fabs? - if unbound then e>=1 -Dan        

    // Mean anamoly
    double M0a = E0 - ecc*sin(E0);  
    if( Etot>0.0 ) M0a = E0 - ecc*sinh(E0); // above unbound

    double Ma = M0a + dt*pow((Ms+Mp)/asmaj, 1.5);//dt*sqrt(mt)/pow(a,1.5);  // advance to new mean anamoly
    if( Etot>0.0 ) Ma = M0a - dt*pow((Ms+Mp)/asmaj, 1.5); //- dt*sqrt(mt)/pow(-a,1.5);

    double E1 = Ma;
    double dfnc = Ma-E1+ecc*sin(E1); //should equal zero M=E-esinE
    // solve for E1 eccentric anamoly
    if( Etot>0.0 ) dfnc = Ma-E1+ecc*sinh(E1);
    while( fabs(dfnc) > 1e-10 ){
      double dfdE = -1.+ecc*cos(E1);
      if( Etot>0.0 ) dfdE = -1.+ecc*cosh(E1);
      double dE = -dfnc/dfdE;
      E1 += dE;
      dfnc = Ma-E1+ecc*sin(E1);
      if( Etot>0.0 ) dfnc = Ma-E1+ecc*sinh(E1);
    }

    rsep = asmaj*(1.-ecc*cos(E1));
    if( Etot>0.0 ) rsep = asmaj*(1.-ecc*cosh(E1));
    double phi1 = 2.*atan2( sqrt(1.+ecc)*sin(E1/2.) , sqrt(1.-ecc)*cos(E1/2.) );
    if( Etot>0.0 )  phi1 = 2.*atan2( sqrt(1.+ecc)*sinh(E1/2.) , sqrt(fabs(1.-ecc))*cosh(E1/2.) );


        
  
    double hh = Ltot*(Ms+Mp)/(Ms*Mp);
    ////double vp0 = hh/(1.+Mp/Ms)/rsep;
    double vp1 = hh/(1.+Ms/Mp)/rsep;
    //double vr = sqrt( -2*Etot/rsep/rsep * ( 0.5*Ltot*Ltot/Etot - 0.5*(Ms+Mp)/Etot * rsep - rsep*rsep  ) );
    double vr = sqrt( fabs( (Ms+Mp)*(2./rsep - 1./asmaj) -hh*hh/rsep/rsep )  ); //- (Ms+Mp)/asmaj/asmaj/asmaj * rsep*rsep); //0 for circ orbit


    //if((Ms+Mp)*(2./rsep - 1./asmaj) - (Ms+Mp)/asmaj/asmaj/asmaj * rsep*rsep < 0.0 ) vr = 0.0;
    if( sin(E1) < 0.0 ) vr = -vr; //???

    thisGravMass->r = rsep/(1.+Ms/Mp);
    ////thisGravMass->omega = pow(rsep,(-3./2.)); 
    //thisGravMass->omega = vp1/rsep/(1.+Ms/Mp);     // vphi/r 
    thisGravMass->vr = vr/(1.+Ms/Mp);
    thisGravMass->omega = hh/rsep/rsep; //fabs(phi0 - phi1)/dt;           //Circ case - Is this a bad idea?
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
  if (sim_GravMassType(theSim)==LIVEBINARY || BINARY){
    //printf("LIVEBINARY");
    double Mp = theGravMasses[0].M; //primary
    double Ms = theGravMasses[1].M; //secondary

    //Here we need to analytically update a, phi, vr from Kepler's equations 
    gravMass_adv_anly(&(theGravMasses[1]), Mp, dt);
    // Or just push the binary together arbitrarily for testing
    //gravMass_adv_arb(theSim, &(theGravMasses[1]), Mp, t, dt); //TESTING 
                                                                       
    theGravMasses[0].r   = theGravMasses[1].r*Ms/Mp;
    theGravMasses[0].omega  = theGravMasses[1].omega; // True for eccentric binary?
    theGravMasses[0].vr  = theGravMasses[1].vr * Ms/Mp;


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
  // IF LIVE BINARY UPDATE MASSES ENERY, ANG MOMENTUM (IS THIS THE BEST PLACE TO DO THIS)?
  if ( (sim_GravMassType(theSim)==LIVEBINARY)  && (time_global > 2.*M_PI*sim_tmig_on(theSim)) ){

    //toal torque gives the torque on the binary due to the gas wrt to r=0 
    double Tot_Trq =  gravMass_total_torque(theGravMasses,0);
    double Tot_Pow =  gravMass_total_power(theGravMasses,0);
    //double Macc0   =  gravMass_Macc(theGravMasses,0);
    //double Macc1   =  gravMass_Macc(theGravMasses,1);
    double Ebin    =  gravMass_E(theGravMasses,1);
    double Lbin    =  gravMass_Ltot(theGravMasses,1); //theGravMasses[1].Ltot;
    //IMPORTANT 0 and 1 Ltot are the same but use 1 only! see adv_anyl() above    

    double dens_scale = ( sim_Mdsk_o_Ms(theSim)/( 1.+1./sim_MassRatio(theSim) ) )/( M_PI*sim_sep0(theSim)*sim_sep0(theSim) );
    if (   time_global <=    2.*M_PI*(sim_tmig_on(theSim) + sim_tramp(theSim))   ){
      dens_scale *= (time_global - 2*M_PI*sim_tmig_on(theSim))/( sim_tramp(theSim)*2*M_PI );
    }      //^REDUNDANT WITH celll-source

       
    //UPDATE E and L
    //ORDER MATTERS FOR UPDATE BELOW - is it correct?
    if (sim_LiveAcc(theSim) == 1){
      double M0   =  gravMass_M(theGravMasses,0);
      double M1   =  gravMass_M(theGravMasses,1);
      double rbh0   =  gravMass_r(theGravMasses,0);
      double rbh1   =  gravMass_r(theGravMasses,1);
      double om     =  gravMass_omega(theGravMasses,1);
      double Mdot0   =  gravMass_Mdot(theGravMasses,0);
      double Mdot1   =  gravMass_Mdot(theGravMasses,1);
      //SCALE accretion in terms of disk density
      double DS_Mdot0 = dens_scale*Mdot0;
      double DS_Mdot1 = dens_scale*Mdot1;
      Lbin  += (DS_Mdot0*rbh0*rbh0* + DS_Mdot1*rbh1*rbh1)*om;
      Ebin  += Ebin*(DS_Mdot0/M0 + DS_Mdot1/M1);
      //UPDATE MASSES
      M0 += DS_Mdot0*dt;
      M1 += DS_Mdot1*dt;
      GravMass_set_M(theGravMasses, 0, M0);
      GravMass_set_M(theGravMasses, 1, M1);
    }
    //else{
    //  // keep track of bin radial velocity for keplerian ecc orbit
    //  double Ebin    =  gravMass_E(theGravMasses,1);
    //  double semaj = 0.5*Ms*Mp/fabs(Ebin);;
    //  double ecc = sim_ecc0(theSim);
    //  double vr = semaj * sqrt(Mtotal/semaj/semaj/semaj) * ecc/sqrt(1.-ecc*ecc) * sin(theGravMasses[1].phi);
    //  theGravMasses[0].vr    = vr/(1.+M0/M1);
    //  theGravMasses[1].vr    = vr/(1.+M1/M0);
    //}


    //due to forces from gas
    Lbin  += Tot_Trq*dt;
    Ebin  += Tot_Pow*dt;

      
    GravMass_set_Ltot(theGravMasses, 1, Lbin);
    GravMass_set_Etot(theGravMasses, 1, Ebin);
   // also update the binary mass and Energy here


  }
}

//PUT THESE IN _set no?
void gravMass_set_Mdot( struct GravMass * theGravMasses, double Mdot, int p){
  theGravMasses[p].Mdot = Mdot;
}

void gravMass_set_Macc( struct GravMass * theGravMasses, double Macc, int p){
  theGravMasses[p].Macc = Macc;
}

void gravMass_set_total_torque( struct GravMass * theGravMasses, double total_torque, int p){
  theGravMasses[p].total_torque = total_torque;
}

void gravMass_set_total_power( struct GravMass * theGravMasses, double total_power, int p){
  theGravMasses[p].total_power = total_power;
}
