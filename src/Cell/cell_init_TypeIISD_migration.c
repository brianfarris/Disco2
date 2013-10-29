#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Sim.h"
#include "../Headers/Face.h"
#include "../Headers/GravMass.h"
#include "../Headers/header.h"

void cell_single_init_TypeIISD_migration(struct Cell *theCell, struct Sim *theSim,int i,int j,int k){

	double rho   = 1.0;
	double Mach   = 10.0;
	
	
	double Gam = sim_GAMMALAW(theSim);
	double fac = 0.0001;

	double Mtotal = 1.0;
	double sep = sim_sep0(theSim);
	double massratio = sim_MassRatio(theSim);
	double M0 = Mtotal/(1.+massratio);
	double M1 = Mtotal/(1.+1./massratio);
	double r0 = M1/Mtotal*sep;
	double r1 = M0/Mtotal*sep;
	
	double eps0 = sim_G_EPS(theSim)*r0;
	double eps1 = sim_G_EPS(theSim)*r1;
	
	double DISK_ALPHA = sim_EXPLICIT_VISCOSITY(theSim);
	
	double rm = sim_FacePos(theSim,i-1,R_DIR);
	double rp = sim_FacePos(theSim,i,R_DIR);
	double r = 0.5*(rm+rp);

			
	double cs;
	cs = pow(sep,(-0.5))/(1.+massratio)/Mach;
	//if (r<1.){
	//	cs = 1./Mach;
	//}else{
	//	cs = sqrt(1./r)/Mach;
	//}
	double xbh0 = r0;
	double xbh1 = -r1;
	double ybh0 = 0.0;
	double ybh1 = 0.0;
			

				
	double phi = theCell->tiph-.5*theCell->dphi;
	double xpos = r*cos(phi);
	double ypos = r*sin(phi);
				
	double dist_bh1 = sqrt((xpos-xbh0)*(xpos-xbh0)+(ypos-ybh0)*(ypos-ybh0));
	double dist_bh2 = sqrt((xpos-xbh1)*(xpos-xbh1)+(ypos-ybh1)*(ypos-ybh1));
				
	double n = sim_PHI_ORDER(theSim);
	double Pot1 = M0/pow( pow(dist_bh1,n) + pow(eps0,n) , 1./n );
	double Pot2 = M1/pow( pow(dist_bh2,n) + pow(eps1,n) , 1./n );
				

	double omega = 1./pow(r,1.5);
	omega *= 1.+3./4.*(sep*sep/r/r*(massratio/((1.+massratio)*(1.+massratio))));
	double O2 = omega*omega; //+ cs*cs/r/r*( 2.*rs*rs/r/r - 3. ); Add Deriv of Pressure term
	double vr;
	omega = sqrt(O2);
	vr = 0.0;

	//if (r<sep){
	//  omega = sqrt(O2);//*pow(r,2.);
	//	vr = 0.0;
	//	//vr = scaled SS;
	//}else{
	//	omega = sqrt(O2);
	//	vr = 0.0;
	//	//vr = scaled SS;
	//}		  
	
	if( rho<sim_RHO_FLOOR(theSim) ) rho = sim_RHO_FLOOR(theSim);
				
	double Pp = cs*cs*rho/Gam;	
	
	
	theCell->prim[RHO] = rho*exp(fac*(Pot1+Pot2)/cs/cs);;
	theCell->prim[PPP] = Pp*exp(fac*(Pot1+Pot2)/cs/cs);
	theCell->prim[URR] = vr;
	theCell->prim[UPP] = omega;
	theCell->prim[UZZ] = 0.0;
	theCell->wiph = 0.0;
	theCell->divB = 0.0;
	theCell->GradPsi[0] = 0.0;
	theCell->GradPsi[1] = 0.0;
	theCell->GradPsi[2] = 0.0;	
}

void cell_init_TypeIISD_migration(struct Cell ***theCells,struct Sim *theSim,struct MPIsetup * theMPIsetup) {
	
  double rho0   = 1.0;
  double Mach   = 10.0;

  double Gam = sim_GAMMALAW(theSim);
  double fac = 0.0001;

  //printf("rbh1: %e, rbh2: %e, mbh1: %e, mbh2: %e\n",gravMass_r(theGravMasses,0),gravMass_r(theGravMasses,1),gravMass_M(theGravMasses,0),gravMass_M(theGravMasses,1));

  double Mtotal = 1.0;
  double sep = sim_sep0(theSim);
  double massratio = sim_MassRatio(theSim);
  double M0 = Mtotal/(1.+massratio);
  double M1 = Mtotal/(1.+1./massratio);
  double r0 = M1/Mtotal*sep;
  double r1 = M0/Mtotal*sep;

  double eps0 = sim_G_EPS(theSim)*r0;
  double eps1 = sim_G_EPS(theSim)*r1;
  
  double DISK_ALPHA = sim_EXPLICIT_VISCOSITY(theSim);

  int i, j,k;
  for (k = 0; k < sim_N(theSim,Z_DIR); k++) {
    for (i = 0; i < sim_N(theSim,R_DIR); i++) {
      double rm = sim_FacePos(theSim,i-1,R_DIR);
      double rp = sim_FacePos(theSim,i,R_DIR);
      double r = 0.5*(rm+rp);

      double cs;
      cs = pow(sep,(-0.5))/(1.+massratio)/Mach;

//      if (r<1.){
 //       cs = 1./Mach;
  //    }else{
   //     cs = sqrt(1./r)/Mach;
    //  }

      double xbh0 = r0;
      double xbh1 = -r1;
      double ybh0 = 0.0;
      double ybh1 = 0.0;

      for (j = 0; j < sim_N_p(theSim,i); j++) {

        double phi = theCells[k][i][j].tiph - 0.5*theCells[k][i][j].dphi;
        double xpos = r*cos(phi);
        double ypos = r*sin(phi);

        double dist_bh1 = sqrt((xpos-xbh0)*(xpos-xbh0)+(ypos-ybh0)*(ypos-ybh0));
        double dist_bh2 = sqrt((xpos-xbh1)*(xpos-xbh1)+(ypos-ybh1)*(ypos-ybh1));

        double n = sim_PHI_ORDER(theSim);
        double Pot1 = M0/pow( pow(dist_bh1,n) + pow(eps0,n) , 1./n );
        double Pot2 = M1/pow( pow(dist_bh2,n) + pow(eps1,n) , 1./n );

	    //---------------------------Shakura Sunyaev---------------------------//
		  //(7.21337*10^23 (1 - 0.0129692 Sqrt[1/R])^(7/10))/R^(3/4)  From Mathematica ".nb" 
		  // For 10^7Msun q=0.01 such that at r_sim=10 Mdisk = Msec 
	//double f = pow(1. - 0.0129692*sqrt(1./r),1./4.);
	double rho = rho0; //pow(r,(-3./4.)) * pow(f,(14./5.));
        double omega = 1./pow(r,1.5);
        omega *= 1.+3./4.*(sep*sep/r/r*(massratio/((1.+massratio)*(1.+massratio))));
	double O2 = omega*omega; //+ cs*cs/r/r*( 2.*rs*rs/r/r - 3. ); Add Deriv of Pressure term
        double vr;
	omega = sqrt(O2);//*pow(r,2.);
	vr = 0.0;

//        if (r<sep){
//	  omega = sqrt(O2);//*pow(r,2.);
//         vr = 0.0;
//          //vr = scaled SS;
//       }else{
//          omega = sqrt(O2);
//          vr = 0.0;
//          //vr = scaled SS;
//        }		  

		//---------------------------Shakura Sunyaev---------------------------//
        /* 
        if (r<3.){
          omega = 0.0;//pow(r,-1.5);
        }
        */
        if( rho<sim_RHO_FLOOR(theSim) ) rho = sim_RHO_FLOOR(theSim);
        
        double Pp = cs*cs*rho/Gam;

        theCells[k][i][j].prim[RHO] = rho*exp(fac*(Pot1+Pot2)/cs/cs);
        theCells[k][i][j].prim[PPP] = Pp*exp(fac*(Pot1+Pot2)/cs/cs);
        theCells[k][i][j].prim[URR] = vr;
        theCells[k][i][j].prim[UPP] = omega;
        theCells[k][i][j].prim[UZZ] = 0.0;
        theCells[k][i][j].wiph = 0.0;
        theCells[k][i][j].divB = 0.0;
        theCells[k][i][j].GradPsi[0] = 0.0;
        theCells[k][i][j].GradPsi[1] = 0.0;
        theCells[k][i][j].GradPsi[2] = 0.0;
      }
    }
  }

}
