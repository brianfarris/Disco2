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
        double massratio = sim_MassRatio(theSim);
        double eps = sim_G_EPS(theSim);// * r0;

	//double OmegaK = pow(r,-1.5);	
	double Omeff = sqrt(pow(dist_bh0*dist_bh0+eps*eps,-1.5)*M0+pow(dist_bh1*dist_bh1+eps*eps,-1.5)*M1);
	double cs = r*Omeff*HoR; //cst MAch
	//double cs = HoR; //cst cs	
	double PoRho = cs*cs/Gam;
	
	
        double rho = c->prim[RHO]; //exponential atmosphere already in here
	//
        c->prim[PPP] = PoRho * rho;
      }
    }
  }
}
