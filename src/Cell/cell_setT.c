#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Sim.h"
#include "../Headers/GravMass.h"
#include "../Headers/header.h"


void cell_setT( struct Cell *** theCells ,struct Sim * theSim, struct GravMass * theGravMasses ){
  double Gam = sim_GAMMALAW(theSim); //DD
  double Mtotal = 1.0;
  //double sep = 1.0;
  double M0 = gravMass_M(theGravMasses,0);
  double M1 = gravMass_M(theGravMasses,1);

  double r0 = gravMass_r(theGravMasses,0);
  double r1 = gravMass_r(theGravMasses,1);

  double phi_bh0 = gravMass_phi(theGravMasses,0);
  double phi_bh1 = gravMass_phi(theGravMasses,1);

  double xbh0 = r0*cos(phi_bh0);
  double ybh0 = r0*sin(phi_bh0);

  double xbh1 = r1*cos(phi_bh1);
  double ybh1 = r1*sin(phi_bh1);
        

  double HoR = sim_HoR(theSim);
  double T0 = 0.02*(HoR/0.1)*(HoR/0.1);
  double r_a = sqrt(1./log(T0*sim_GAMMALAW(theSim)/(HoR*HoR)));
  
  int i,j,k;
  for( k=0 ; k<sim_N(theSim,Z_DIR) ; ++k ){
    double zm = sim_FacePos(theSim,k-1,Z_DIR);
    double zp = sim_FacePos(theSim,k,Z_DIR);
    for( i=0 ; i<sim_N(theSim,R_DIR) ; ++i ){
      double rm = sim_FacePos(theSim,i-1,R_DIR);
      double rp = sim_FacePos(theSim,i,R_DIR);
      double r = 0.5*(rm+rp);
      for( j=0 ; j<sim_N_p(theSim,i) ; ++j ){
        struct Cell * c = &(theCells[k][i][j]);

	
        double phi = c->tiph - 0.5*c->dphi;
        double xpos = r*cos(phi);
        double ypos = r*sin(phi);

        //printf("xpos: %e, ypos: %e\n",xpos,ypos);

        double dist_bh0 = sqrt((xpos-xbh0)*(xpos-xbh0)+(ypos-ybh0)*(ypos-ybh0));
        double dist_bh1 = sqrt((xpos-xbh1)*(xpos-xbh1)+(ypos-ybh1)*(ypos-ybh1));
        //printf("dist_bh0: %e, dist_bh1: %e\n",dist_bh0,dist_bh1);


        double Mtotal = 1.0;
        //double sep = sim_sep0(theSim);
        double massratio = sim_MassRatio(theSim);
        //double M0 = Mtotal/(1.+massratio);
        //double M1 = Mtotal/(1.+1./massratio);
        //double r0 = sep*M1/Mtotal;
        //double r1 = sep*M0/Mtotal;
        double eps0 = sim_G_EPS(theSim) * r0;
        double eps1 = sim_G_EPS(theSim) * r1;

       
        //printf("dist_bh0: %e, dist_bh1: %e, xbh0: %e, ybh0: %e, xbh1: %e, ybh1: %e\n",dist_bh0,dist_bh1,xbh0,ybh0,xbh1,ybh1);
	//double cs0,cs1;
        //double PoRho0,PoRho1,PoRho;
        /*
	double vp0,vp1;
	// before the cutoffs around each BH were set to 0.05, but in milos ICs were set to 0.5? I set them to rsink here
        if (dist_bh0>sim_Rsink0(theSim)){
		vp0 = pow(dist_bh0,-0.5) * sqrt(M0) ;
        	cs0 = HoR*vp0;
                PoRho0 = cs0*cs0/sim_GAMMALAW(theSim);
        } else{
	        vp0 = pow(sim_Rsink0(theSim),-0.5);
		cs0 = HoR*vp0;
                PoRho0 = cs0*cs0/sim_GAMMALAW(theSim);
                //PoRho = T0 * exp(-r*r/(r_a*r_a)); 
	}
        if (dist_bh1>sim_Rsink1(theSim)){
		vp1 = pow(dist_bh1,-0.5) * sqrt(M1) ;
        	cs1 = HoR*vp1;
                PoRho1 = cs1*cs1/sim_GAMMALAW(theSim);
        } else{
	        vp1 = pow(sim_Rsink1(theSim),-0.5);
		cs1 = HoR*vp1;
                PoRho1 = cs1*cs1/sim_GAMMALAW(theSim);
                //PoRho = T0 * exp(-r*r/(r_a*r_a)); 
	}


	

        PoRho = PoRho0+PoRho1;
	*/





	


        //cs0 = pow((dist_bh0*dist_bh0 + 0.05*0.05),-0.5) * sqrt(M0) * HoR;
        //cs1 = pow((dist_bh1*dist_bh1 + 0.05*0.05),-0.5) * sqrt(M1) * HoR;
        //PoRho0 =  cs0*cs0/Gam;
        //PoRho1 =  cs1*cs1/Gam;
	//double cs2 = (cs0*cs0 + cs1*cs1);

	//double fac = 0.0001;
	//double n = sim_PHI_ORDER(theSim);
	//double Pot0 = M0/pow( pow(dist_bh0,n) + pow(eps0,n) , 1./n );
	//double Pot1 = M1/pow( pow(dist_bh1,n) + pow(eps1,n) , 1./n );
	

	double OmegaK = pow(r,-1.5);
	double cs = r*OmegaK*HoR; //cst MAch
	//double cs = HoR; //cst cs
	
	double PoRho = cs*cs/Gam;//*exp(fac*(Pot0+Pot1)/cs/cs);

	//PoRho = (PoRho0+PoRho1)*exp(fac*(Pot0+Pot1)/cs2);
	
	
        double rho = c->prim[RHO];
        c->prim[PPP] = PoRho * rho;
      }
    }
  }
}
