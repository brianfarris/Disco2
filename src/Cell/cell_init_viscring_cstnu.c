#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Sim.h"
#include "../Headers/Face.h"
#include "../Headers/GravMass.h"
#include "../Headers/header.h"

void cell_single_init_viscring_cstnu(struct Cell *theCell, struct Sim *theSim,int i,int j,int k){
	double R0     = 1.0;
	double GM     = 1.0;
	double rho0   = 1.0;
	double Mach   = 100.0;
	
	
	

	double nu = sim_EXPLICIT_VISCOSITY(theSim); //DISK_ALPHA*(1.0/Mach)*(1.0/Mach);
	if (nu<=0.){
		nu=0.;
	}
	
	//double tau = 12.*nu*t/R0^2;
	double tau;
	double tau0 = 0.032;
	if (fabs(nu)<0.000001) {
       	tau = tau0;
	}else {
		double t    = tau0*R0*R0/(12.*nu) + time_global;
		double tau  = 12.*nu/(R0*R0) * t;
	}
		
	
	double Gam = sim_GAMMALAW(theSim);

	double rm = sim_FacePos(theSim,i-1,R_DIR);
	double rp = sim_FacePos(theSim,i,R_DIR);
	double r = 0.5*(rm+rp);
	
	double cs;
	//if (r<1.){
	cs = 1./Mach;
	// }else{
	//cs = sqrt(1./r)/Mach;
	//}
	
	double dx = r/R0;
	double argM = 2.*dx/tau;
	
	double I14  = exp(argM)/sqrt(2.*M_PI*argM)*( 1.+3./(32.*argM) );
	
	double In34 = exp(argM)/sqrt(2.*M_PI*argM)*( 1.-5./(32.*argM) );
	
	double I54  = exp(argM)/sqrt(2.*M_PI*argM)*( 1.-21./(32.*argM) );
	
	double rho = rho0* 1./tau * pow(dx,-0.25) *exp( -(1.+dx*dx)/tau )*I14;
	double vp  = sqrt(GM/r);
	double vr  = -3.*nu/R0 * (  2./tau*(In34/I14 - r/R0)   );
	double omega = vp/r;
	
	//double n = sim_PHI_ORDER(theSim);
	//double Pot = 1./pow( pow(r,n) + pow(eps,n) , 1./n );
	
	if( rho<sim_RHO_FLOOR(theSim) ) rho = sim_RHO_FLOOR(theSim);
	double Pp  = rho*cs*cs/Gam;
	//double Pp = 1.0;
	
	theCell->prim[RHO] = rho;
	theCell->prim[PPP] = Pp;
	theCell->prim[URR] = vr;
	theCell->prim[UPP] = omega;
	theCell->prim[UZZ] = 0.0;
	theCell->wiph = 0.0;
		
}



void cell_init_viscring_cstnu(struct Cell ***theCells,struct Sim *theSim,struct MPIsetup * theMPIsetup) {
  double R0     = 1.0;
  double GM     = 1.0;
  double rho0   = 1.0;
  double Mach   = 100.0;


  double nu = sim_EXPLICIT_VISCOSITY(theSim); //DISK_ALPHA*(1.0/Mach)*(1.0/Mach);
  if (nu<=0.){
	nu=0.;
  }
  //double tau = 12.*nu*t/R0^2;
  double tau = 0.032;

  double Gam = sim_GAMMALAW(theSim);
  //double r_planet = 0.1;
  //double eps = sim_G_EPS(theSim)*r_planet;
  //double cp = 1./Mach;
  //double fac = 0.0001;



  int i,j,k;
  for (k = 0; k < sim_N(theSim,Z_DIR); k++) {
    for (i = 0; i < sim_N(theSim,R_DIR); i++) {
      double rm = sim_FacePos(theSim,i-1,R_DIR);
      double rp = sim_FacePos(theSim,i,R_DIR);
      double r = 0.5*(rm+rp);
		
      double cs;
      //if (r<1.){
	  cs = 1./Mach;
     // }else{
	  //cs = sqrt(1./r)/Mach;
      //}
				
	  double dx = r/R0;
	  double argM = 2.*dx/tau;
		
	  double I14  = exp(argM)/sqrt(2.*M_PI*argM)*( 1.+3./(32.*argM) );
		
	  double In34 = exp(argM)/sqrt(2.*M_PI*argM)*( 1.-5./(32.*argM) );
		
	  double I54  = exp(argM)/sqrt(2.*M_PI*argM)*( 1.-21./(32.*argM) );
		
	  double rho = rho0 * 1./tau * pow(dx,-0.25) *exp( -(1.+dx*dx)/tau )*I14;
	  double vp  = sqrt(GM/r);
	  double vr  = -3.*nu/R0 * (  2./tau*(In34/I14 - r/R0)   );
	  double omega = vp/r;
		
	  //double n = sim_PHI_ORDER(theSim);
	  //double Pot = 1./pow( pow(r,n) + pow(eps,n) , 1./n );
		
	  if( rho<sim_RHO_FLOOR(theSim) ) rho = sim_RHO_FLOOR(theSim);
	  double Pp  = rho*cs*cs/Gam;
	  //double Pp = 1.0;

      for (j = 0; j < sim_N_p(theSim,i); j++) {
		  theCells[k][i][j].prim[RHO] = rho;//*exp(fac*(Pot)/cp/cp);
		  theCells[k][i][j].prim[PPP] = Pp;//*exp(fac*(Pot)/cp/cp);
          theCells[k][i][j].prim[URR] = vr;
          theCells[k][i][j].prim[UPP] = omega;
          theCells[k][i][j].prim[UZZ] = 0.0;
          theCells[k][i][j].wiph = 0.0;
      }
    }
  }

}
