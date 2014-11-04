#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Sim.h"
#include "../Headers/Face.h"
#include "../Headers/GravMass.h"
#include "../Headers/header.h"

void cell_single_init_SimpleKepler(struct Cell *theCell, struct Sim *theSim,int i,int j,int k){

	double rho   = 1.0;
	double Gam = sim_GAMMALAW(theSim);
	
	double sep = sim_sep0(theSim);
	double HoR = sim_HoR(theSim);
	double Mach   = 1./HoR;                   //For uniform Mach	


	double rm = sim_FacePos(theSim,i-1,R_DIR);
	double rp = sim_FacePos(theSim,i,R_DIR);
	double r = 0.5*(rm+rp);
	


	double OmegaK = pow(r,-1.5);
	double omega = OmegaK;

	//double cs = r*OmegaK/Mach;
	double cs = 1./Mach;


	double PoRho = cs*cs/Gam;   //for any gamma law eqn of state (adiabatic or iso)				


        double disk_nu;
	double DISK_ALPHA = sim_EXPLICIT_VISCOSITY(theSim);
        if (sim_VISC_CONST(theSim)==1){
          disk_nu = DISK_ALPHA;
	}else{
	  printf("THESE ICS ARE NOT CONSISTENT WITH NON-CONST VISCOSITY");
        }

   

	double vr;
        if (DISK_ALPHA > 0.0){
          vr = -3./2.*disk_nu/r;
	    }else{
          vr = 0.0;
        }


	
	if( rho<sim_RHO_FLOOR(theSim) ) rho = sim_RHO_FLOOR(theSim);
				

	double Pp = PoRho * rho;


	double passive_scalar;
	if (r>1.3 && r<1.5){
	  passive_scalar = 1.;	  
	}else{
	  passive_scalar = 0.000000001;
	}


	theCell->prim[RHO] = rho;
	theCell->prim[PPP] = Pp;
	theCell->prim[URR] = vr;
	theCell->prim[UPP] = omega - sim_rOm_a(theSim,r,sep)/r; //initial sep set to 1.
	theCell->prim[UZZ] = 0.0;
	theCell->wiph = 0.0;//r*OmegaK;
>>>>>>> fdcf1c6801e932c8fa1f2a128c51a205b4fc4d44
	if(sim_NUM_C(theSim)<sim_NUM_Q(theSim)) theCell->prim[sim_NUM_C(theSim)] = passive_scalar;	
}

void cell_init_SimpleKepler(struct Cell ***theCells,struct Sim *theSim,struct MPIsetup * theMPIsetup) {
	
  double rho0   = 1.0;

  double Gam = sim_GAMMALAW(theSim);

  //printf("rbh1: %e, rbh2: %e, mbh1: %e, mbh2: %e\n",gravMass_r(theGravMasses,0),gravMass_r(theGravMasses,1),gravMass_M(theGravMasses,0),gravMass_M(theGravMasses,1));

  double sep = sim_sep0(theSim);
  double HoR = sim_HoR(theSim);   // for uniform Mach
  double Mach   = 1./HoR;
  
  double DISK_ALPHA = sim_EXPLICIT_VISCOSITY(theSim);

  int i, j,k;
  for (k = 0; k < sim_N(theSim,Z_DIR); k++) {
    for (i = 0; i < sim_N(theSim,R_DIR); i++) {
      double rm = sim_FacePos(theSim,i-1,R_DIR);
      double rp = sim_FacePos(theSim,i,R_DIR);
      double r = 0.5*(rm+rp);


      for (j = 0; j < sim_N_p(theSim,i); j++) {

        double phi = theCells[k][i][j].tiph - 0.5*theCells[k][i][j].dphi;
        double xpos = r*cos(phi);
        double ypos = r*sin(phi);


	double OmegaK = 1./pow(r,1.5);

	//double cs = r*OmegaK/Mach;
	double cs = 1./Mach;

	double PoRho = cs*cs/Gam;

	double rho = rho0;



        double disk_nu;
        double DISK_ALPHA = sim_EXPLICIT_VISCOSITY(theSim);
        if (sim_VISC_CONST(theSim)==1){
          disk_nu = DISK_ALPHA;
        }else{
	  printf("THESE ICS ARE NOT CONSISTENT WITH NON-CONST VISCOSITY");
        }




        double omega = OmegaK;


        double vr;
	if (DISK_ALPHA > 0.0){
          vr = -3./2.*disk_nu/r; 
        }else{
          vr = 0.0;
        }
	


        if( rho<sim_RHO_FLOOR(theSim) ) rho = sim_RHO_FLOOR(theSim);
        

	double Pp = PoRho * rho;



	double passive_scalar;
        if (r>1.3 && r<1.5){
          passive_scalar = 1.;
	}else{
	  passive_scalar = 0.000000001;;
	}





        theCells[k][i][j].prim[RHO] = rho;
        theCells[k][i][j].prim[PPP] = Pp;
        theCells[k][i][j].prim[URR] = vr;
        theCells[k][i][j].prim[UPP] = omega -  sim_rOm_a(theSim,r,sep)/r; //initial sep set to 1. 
        theCells[k][i][j].prim[UZZ] = 0.0;
        theCells[k][i][j].wiph = 0.0; //r*OmegaK;
	if(sim_NUM_C(theSim)<sim_NUM_Q(theSim)) theCells[k][i][j].prim[sim_NUM_C(theSim)] = passive_scalar;
      }
    }
  }

}
