#define PLANET_PRIVATE_DEFS
#include <stdlib.h>
#include <math.h>
#include "../Headers/GravMass.h"
#include "../Headers/Sim.h"
#include "../Headers/Cell.h"
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
	theGravMasses[p].RK_vr  = theGravMasses[p].vr;
  }
}


/* PAUL's DISCO Implementation adjusted for DISCO2 (DO NOT USE YET)
void gravMass_adv_anly( struct GravMass * thisGravMass , double MM , double dt ){
	
	double m = thisGravMass->M;      //The mass of the secondary
	double mt = MM/pow(1.+m/MM,2.);  // Mprim/(1+q)^2 
	double E = thisGravMass->E/m/(1.+m/MM); // Etot/(Mbin)
	double L = thisGravMass->L/m;           // Lsec/Msec
	
	double r = thisGravMass->r;
	double phi = thisGravMass->phi;
	double vr = thisGravMass->vr;

	double a = -.5*mt/E;               //This is -MsMp/(2E*(1+q)) but shouldnt it be a = -MsMp/(2E)? -Dan
	double e = sqrt( 1. + 2.*L*L*E/mt/mt ); // This is sqrt(1+2L^2E/(MsMp)^2*(1+q)) but shouldnt the (q+1 not be there?)
	//if( 1. + 2.*L*L*E < 0.0 ) e = 0.0;   //missing mt^2???
	if( 1. + 2.*L*L*E/mt/mt < 0.0 ) e = 0.0;
	
	double cosE0 = (1.-r/a)/e;    // from r=a(1-ecos(E))
	if( e <= 0.0 ) cosE0 = 0.0;   // circular orbit or imaginary e -> something weird?
	double sinE0 = sqrt(fabs(1.-cosE0*cosE0));
	if( vr < 0.0 ) sinE0 = -sinE0; //???
	
	double E0 = atan2( sinE0 , cosE0 );
	if( E>0.0 ) E0 = log( sinE0 + cosE0 );
	
	double phi0 = 2.*atan2( sqrt(1.+e)*sin(E0/2.) , sqrt(1.-e)*cos(E0/2.) );   //standard Newtonian formulae for phi
	//if( E>0.0 ) phi0 = 2.*atan2( sqrt(1.+e)*sinh(E0/2.) , sqrt(fabs(1.-e))*cosh(E0/2.) ); //Paul
	if( E>0.0 ) phi0 = 2.*atan2( sqrt(e+1.)*sinh(E0/2.) , sqrt(e-1.)*cosh(E0/2.) );//above unbound; why use fabs? - if unbound then e>=1 -Dan 
	
	double M0 = E0 - e*sin(E0);   // Newtonain formulae for t(=M here) in terms of Eccentric anomaly E and eccentricity e
	if( E>0.0 ) M0 = E0 - e*sinh(E0); // above unbound
	
	double M = M0 + dt*sqrt(mt)/pow(a,1.5);  // advance t but why with Omega = sqrt(mt/a^3) ??? 
	if( E>0.0 ) M = M0 - dt*sqrt(mt)/pow(-a,1.5);
	
	double E1 = M;
	double df = M-E1+e*sin(E1);  // ???
	if( E>0.0 ) df = M-E1+e*sinh(E1);
	while( fabs(df) > 1e-12 ){
		double dfdE = -1.+e*cos(E1);
		if( E>0.0 ) dfdE = -1.+e*cosh(E1);
		double dE = -df/dfdE;
		E1 += dE;
		df = M-E1+e*sin(E1);
		if( E>0.0 ) df = M-E1+e*sinh(E1);
	}
	
	r = a*(1.-e*cos(E1));    // now set teh new r and phi
	if( E>0.0 ) r = a*(1.-e*cosh(E1));
	double phi1 = 2.*atan2( sqrt(1.+e)*sin(E1/2.) , sqrt(1.-e)*cos(E1/2.) );
	if( E>0.0 )  phi1 = 2.*atan2( sqrt(1.+e)*sinh(E1/2.) , sqrt(fabs(1.-e))*cosh(E1/2.) );
	
	//*rt = r;
	//*phit = phi + (phi1-phi0);
	
	thisGravMass->r = r;
	thisGravMass->phi = phi + (phi1-phi0);
	
	double vr2 = sqrt( 2.*E + 2./r - L*L/r/r ); //???
	
	if( 2.*E + 2./r - L*L/r/r < 0.0 ) vr2 = 0.0;
	if( sin(E1) < 0.0 ) vr2 = -vr2; //???
	
	//*vrt = vr2;
	thisGravMass->vr = vr2;
	
}
*/ 
///////////////////////////////////////

 
void gravMass_adv_anly( struct GravMass * thisGravMass , double Mp , double dt ){
	
	double Ms = thisGravMass->M;  //The mass of the secondary
	double mt = Mp/pow(1.+Ms/Mp,2.);            // ??
	double E = thisGravMass->E;   // Etot
	double Ls = thisGravMass->L;           // Lsec/Msec
	//double Lp = Ls*Ms/Mp; // Omega the same right? 3 1/q's to convert 1M and 2 r's works out to one Ms/Mp
	double Ltot = Ls*(1.+Ms/Mp);
	
	double r1 = thisGravMass->r;
	double phi = thisGravMass->phi;
	double vr1 = thisGravMass->vr;
	
	double rt = r1 + r1*Ms/Mp;
	
	//either update with Delta L 
	double a = Ltot*Ltot/(Ms*Ms*Mp*Mp)*(Ms+Mp);
	
	// Or update with Delta E
	//double a = -.5*Ms*Mp/Enew;               // E = -G(M0+M1)/(2a) for elliptical orbit

	
	// for e=0 only
	double e=0.0;
	double cosE0 = 0.0;   // circular orbit or imaginary e -> something weird?
	double sinE0 = sqrt(fabs(1.-cosE0*cosE0));
	if( vr1 < 0.0 ) sinE0 = -sinE0; //???
	
	double E0 = atan2( sinE0 , cosE0 );
	if( E>0.0 ) E0 = log( sinE0 + cosE0 );
	
	double phi0 = 2.*atan2( sqrt(1.+e)*sin(E0/2.) , sqrt(1.-e)*cos(E0/2.) );   //standard Newtonian formulae for phi
	//if( E>0.0 ) phi0 = 2.*atan2( sqrt(1.+e)*sinh(E0/2.) , sqrt(fabs(1.-e))*cosh(E0/2.) ); //Paul
	if( E>0.0 ) phi0 = 2.*atan2( sqrt(e+1.)*sinh(E0/2.) , sqrt(e-1.)*cosh(E0/2.) );//above unbound; why use fabs? - if unbound then e>=1 -Dan 
	
	double M0 = E0 - e*sin(E0);   // Newtonian formulae for t(=M here) in terms of Eccentric anomaly E and eccentricity e
	if( E>0.0 ) M0 = E0 - e*sinh(E0); // above unbound
	
	double M = M0 + dt*sqrt(mt)/pow(a,1.5);  // advance t 
	if( E>0.0 ) M = M0 - dt*sqrt(mt)/pow(-a,1.5);
	
	double E1 = M;
	double df = M-E1+e*sin(E1);
	if( E>0.0 ) df = M-E1+e*sinh(E1);
	while( fabs(df) > 1e-12 ){
		double dfdE = -1.+e*cos(E1);
		if( E>0.0 ) dfdE = -1.+e*cosh(E1);
		double dE = -df/dfdE;
		E1 += dE;
		df = M-E1+e*sin(E1);
		if( E>0.0 ) df = M-E1+e*sinh(E1);
	}
	
	rt = a*(1.-e*cos(E1));
	if( E>0.0 ) rt = a*(1.-e*cosh(E1));
	double phi1 = 2.*atan2( sqrt(1.+e)*sin(E1/2.) , sqrt(1.-e)*cos(E1/2.) );
	if( E>0.0 )  phi1 = 2.*atan2( sqrt(1.+e)*sinh(E1/2.) , sqrt(fabs(1.-e))*cosh(E1/2.) );
	
	
	thisGravMass->r = rt/(1.+Ms/Mp);
	double vr = sqrt((Ms+Mp)*(2./rt - 1./a) - (Ms+Mp)/a/a/a * rt*rt); //0 for circ orbit

	if((Ms+Mp)*(2./rt - 1./a) - (Ms+Mp)/a/a/a * rt*rt < 0.0 ) vr = 0.0;
	if( sin(E1) < 0.0 ) vr = -vr; //???
	
	//*vrt = vr2;
	thisGravMass->vr = vr/(1.+Ms/Mp);
	thisGravMass->omega = pow(a,(-3./2.)); 		 //Is this a bad idea?
	
}






//Analytic Kepler
void gravMass_move( struct Sim * theSim, struct GravMass * theGravMasses, double dt ){
	if (sim_GravMassType(theSim)==LIVEBINARY){
		//printf("LIVEBINARY");
		double Mp = theGravMasses[0].M;
		double Ms = theGravMasses[1].M;

	
		//Here we need to analytically update a, phi, vr from Kepler's equations
		gravMass_adv_anly(&(theGravMasses[1]), Mp, dt);
	
		//theGravMasses[0].phi = theGravMasses[1].phi + M_PI;
		theGravMasses[0].r   = theGravMasses[1].r*Ms/Mp;
	    theGravMasses[0].vr  = theGravMasses[1].vr*Ms/Mp;
		theGravMasses[0].omega  = theGravMasses[1].omega;
		
		theGravMasses[1].phi += theGravMasses[1].omega*dt; //Instead of updating in adv_anly  
		theGravMasses[0].phi = theGravMasses[1].phi + M_PI;
	}else{
		theGravMasses[0].phi += theGravMasses[0].omega*dt;
		theGravMasses[1].phi += theGravMasses[1].omega*dt;
	}
}
 



//If Binary with CoM off from origin of coordinates - then we need these angles
/*
void gravMass_beta(struct Sim * theSim, struct GravMass * theGravMasses, double * betap, double * betas){
	double q    = sim_MassRatio(theSim);
	double phip = theGravMasses[0].phi;
	double phis = theGravMasses[1].phi;
	
	*betas = (1.+1./q)/(2. + q + 1./q) * (phis-phip);
	*betap = (phis-phip)*(1. - ( (1.+1./q)/(2. + q + 1./q)) );
}

void gravMass_gam(struct Sim * theSim, struct GravMass * theGravMasses, double * gamp, double * gams){
	double q    = sim_MassRatio(theSim);
	double rp   = theGravMasses[0].r;
	double rs   = theGravMasses[1].r;
	double phip   = theGravMasses[0].phi;
	double phis   = theGravMasses[1].phi;
	double abin = sqrt(rs*rs + rp*rp - 2.*rs*rp*cos(phis-phip));

	double betap, betas;
	gravMass_beta(theSim, theGravMasses, &betap, &betas);
	
	//check + or - below and when does gam need to change sign?
	*gamp = acos( (-sin(betap)*sin(betap)  + cos(betap) * sqrt(  abin/(1.+1./q)*abin/(1.+1./q) - rp*rp*sin(betap)*sin(betap) ))/(abin/(1.+1./q)) );
	*gams = acos( (-sin(betas)*sin(betas)  + cos(betas) * sqrt(  abin/(1.+q)*abin/(1.+q) - rs*rs*sin(betas)*sin(betas) ))/(abin/(1.+q)) );
}
 */

/*
//Integrate Orbits Full Live Bin
void gravMass_move( struct Sim * theSim, struct GravMass * theGravMasses, double dt , double RK){
	if (sim_GravMassType(theSim)==LIVEBINARY){
		//printf("LIVEBINARY");
		int i;
		for( i=0 ; i<sim_NumGravMass(theSim) ; ++i ){
			
			if (RK==0.0){//The first step of RK, RK=0.0, dt=dt
				//Full Live Bin
				theGravMasses[i].r   += theGravMasses[i].vr*dt;
				theGravMasses[i].phi += theGravMasses[i].omega*dt; 
			}
			
	
			if (RK==0.5){//The second step of RK, RK=0.5, dt=dt/2
				double rhlf = theGravMasses[i].r;  // so that omega is updated at the correct r
				double Fr = theGravMasses[i].Fr;
				double Fp = theGravMasses[i].Fp; 
				double M = theGravMasses[i].M;
				//theGravMasses[i].r   += theGravMasses[i].vr*dt;   // advance r and phi to n+1 values  
				//theGravMasses[i].phi += theGravMasses[i].omega*dt; 
				theGravMasses[i].r   = theGravMasses[i].RK_r   + theGravMasses[i].vr*dt*2.0;
				theGravMasses[i].phi = theGravMasses[i].RK_phi + theGravMasses[i].omega*dt*2.0; 
			    theGravMasses[i].vr     += (Fr/M)*dt;
			    theGravMasses[i].omega  += (Fp/(M*rhlf))*dt;  // advance velocities to N=1 values after update at midpoint
			}
			
			//theGravMasses[i].r   = (theGravMasses[i].r   + theGravMasses[i].vr*dt )   * (2.0*(RK-0.5)) + (theGravMasses[i].RK_r   + theGravMasses[i].vr*dt )   * (4.0*RK);
			//theGravMasses[i].phi = (theGravMasses[i].phi + theGravMasses[i].omega*dt) * (2.0*(RK-0.5)) + (theGravMasses[i].RK_phi + theGravMasses[i].omega*dt )* (4.0*RK);

			
			//theGravMasses[i].r   += theGravMasses[i].vr*dt;
			//theGravMasses[i].phi += theGravMasses[i].omega*dt;
			
		}
	}else{
		theGravMasses[0].phi += theGravMasses[0].omega*dt;
		theGravMasses[1].phi += theGravMasses[1].omega*dt;
	}
}
////////////////////////////////////////////////////////////////////
*/





void gravMass_update_RK( struct Cell *** theCells, struct GravMass * theGravMasses,struct Sim * theSim, double RK, double dt){
  int i;
  for( i=0 ; i<sim_NumGravMass(theSim) ; ++i ){
    while( theGravMasses[i].phi-theGravMasses[i].RK_phi >  M_PI ) theGravMasses[i].RK_phi += 2.*M_PI;
    while( theGravMasses[i].phi-theGravMasses[i].RK_phi < -M_PI ) theGravMasses[i].RK_phi -= 2.*M_PI;
    theGravMasses[i].r     = (1.0-RK)*theGravMasses[i].r   + RK*theGravMasses[i].RK_r;
    theGravMasses[i].phi   = (1.0-RK)*theGravMasses[i].phi + RK*theGravMasses[i].RK_phi;
    theGravMasses[i].M     = (1.0-RK)*theGravMasses[i].M   + RK*theGravMasses[i].RK_M;
	theGravMasses[i].omega = (1.0-RK)*theGravMasses[i].omega  + RK*theGravMasses[i].RK_omega;
	theGravMasses[i].L     = (1.0-RK)*theGravMasses[i].L   + RK*theGravMasses[i].RK_L;
	theGravMasses[i].E     = (1.0-RK)*theGravMasses[i].E   + RK*theGravMasses[i].RK_E;
	theGravMasses[i].vr    = (1.0-RK)*theGravMasses[i].vr  + RK*theGravMasses[i].RK_vr;
	if (sim_GravMassType(theSim)==LIVEBINARY){
		// calculate forces and velocities at given position and time in RK timestep  (for direct orbit integration)
		//if (RK==0.5){  
		//	cell_gravMassForcePlanets( theSim, theCells, theGravMasses ); // calculate forces on each planet from back reaction of disk (gravMass_update_RK() uses these)
		//}
		cell_gravMassForcePlanets( theSim, theCells, theGravMasses );
		
		//double M  = theGravMasses[i].M;
		double r  = theGravMasses[i].r;
		double vr = theGravMasses[i].vr;
		//double vp = theGravMasses[i].L/M/r;
		double vp = theGravMasses[i].omega*r;
				
		double Fr = 0.0;//theGravMasses[i].Fr; 
		double Fp = theGravMasses[i].Fp;
		
		
		theGravMasses[i].L      += r*Fp*dt;// * (4.0*RK);            // first step RK=0.0 and no update (for direct orbit integration)
		//Keep track of a total E use 0 arbitrarily
		theGravMasses[0].E      += (Fp*vp + Fr*vr)*dt;// * (4.0*RK); // on second step RK=0.5 and the velocities will be the midpoint value (for direct orbit integration)
		theGravMasses[1].E      += (Fp*vp + Fr*vr)*dt; // Keep E as the total ENERGY
		//theGravMasses[i].vr     += (Fr/M)*dt * (2.0*RK);
		//theGravMasses[i].omega  += (Fp/(r*M))*dt * (2.0*RK);
    
    }

  }
}
