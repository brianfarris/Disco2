#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Sim.h"
#include "../Headers/Face.h"
#include "../Headers/GravMass.h"
#include "../Headers/header.h"

void cell_single_init_MMRingCav(struct Cell *theCell, struct Sim *theSim,int i,int j,int k){


        // MM params
        double rho0   = 1.0;
	double rs = 10.;
	double delta_exp   = 3.0;
	double xi_exp = 2.0;
	
	double Gam = sim_GAMMALAW(theSim);
	double fac = 0.0001;

	double Mtotal = 1.0;
	double sep = sim_sep0(theSim);
	double massratio = sim_MassRatio(theSim);
	double M0 = Mtotal/(1.+massratio);
	double M1 = Mtotal/(1.+1./massratio);
	double r0 = sep*M1/Mtotal;
	double r1 = sep*M0/Mtotal;
	double phi0 = 0.0;
	double phi1 = M_PI; 
	double eps0 = sim_G_EPS(theSim) * r0;
        double eps1 = sim_G_EPS(theSim) * r1;	
	double eps = eps0/r0;

	if (massratio == 0.0){
	  sep = 0.0;
	  M1 = 0.0;
	  r0 = 0.0;
	  r1 = 0.0;	 
	  eps0 = sim_G_EPS(theSim);
	  eps1 = eps0;
	}


	double HoR = sim_HoR(theSim);
	double Mach   = sqrt(Gam)*1./HoR;                   //For uniform Mach	


	double rm = sim_FacePos(theSim,i-1,R_DIR);
	double rp = sim_FacePos(theSim,i,R_DIR);
	double r = 0.5*(rm+rp);

	
	double xbh0 = r0;
	double xbh1 = -r1;
	double ybh0 = 0.0;
	double ybh1 = 0.0;
			
				
	double phi = theCell->tiph-.5*theCell->dphi;
	double xpos = r*cos(phi);
	double ypos = r*sin(phi);
				
	double dist_bh0 = sqrt((xpos-xbh0)*(xpos-xbh0)+(ypos-ybh0)*(ypos-ybh0));
	double dist_bh1 = sqrt((xpos-xbh1)*(xpos-xbh1)+(ypos-ybh1)*(ypos-ybh1));
	
	//-------------------------
	//SET CS
	//------------------------
	double OmegaK = pow(r,-1.5);
	double Omeff = sqrt(pow(dist_bh0*dist_bh0+eps*eps,-1.5)*M0+pow(dist_bh1*dist_bh1+eps*eps,-1.5)*M1); // For HS EQL in binary potential
	double Veff = sqrt(pow(dist_bh0*dist_bh0+eps*eps,-0.5)*M0+pow(dist_bh1*dist_bh1+eps*eps,-0.5)*M1); //Virial Thm

	//double cs = r*Omeff/Mach;
	double cs = Veff/Mach;
	double cs2 = cs*cs;
	double PoRho = cs2/Gam;   //for any gamma law eqn of state (adiabatic or iso)

				
	double n = sim_PHI_ORDER(theSim);
	double Pot0 = M0/pow( pow(dist_bh0,n) + pow(eps0,n) , 1./n );
	double Pot1 = M1/pow( pow(dist_bh1,n) + pow(eps1,n) , 1./n );
				
	//-------------------------
	//SET dnu/dr to compute VR
	//------------------------
	double drcs2 = -1./(r*r)/Mach/Mach; //WRONG FOR Omeff
        double disk_nu;
	double dr_nu;
	double DISK_ALPHA = sim_EXPLICIT_VISCOSITY(theSim);
	double dr_Omeff;
	double numer;
        if (sim_VISC_CONST(theSim)==1){
          disk_nu = DISK_ALPHA;
	  dr_nu   = 0.0;
	}else{
          disk_nu = DISK_ALPHA*cs2/Omeff;
	  //
	  numer = ( (r-r0 * cos(phi-phi0)) * pow(dist_bh0,-5.) )/(1.+massratio) + (  (r-r1 * cos(phi-phi1)) * pow(dist_bh1,-5.) )/(1.+1./massratio);
	  //
	  dr_Omeff = - 1.5* numer / sqrt( pow(dist_bh0,-3.)/(1.+massratio) + pow(dist_bh1,-3.)/(1.+1./massratio)   );
	  //
	  dr_nu = DISK_ALPHA * ( drcs2/Omeff - cs2/Omeff/Omeff * dr_Omeff);
        }

	//-------------------------
	//SET Density, Omega
	//------------------------
	double rho = rho0*( pow(rs/r,delta_exp) )*exp(-pow(rs/r,xi_exp) );
        double omega = sqrt(1.)/pow(r,1.5);
        omega *= ( 1.+3./4.*(sep*sep/r/r*(massratio/((1.+massratio)*(1.+massratio)))) );
	double O2 = omega*omega + cs*cs/r/r*( 2.*rs*rs/r/r - 3. );
        double vr;
        if (r<3.0){
          omega = 1./pow(r,1.5);//sqrt(O2);
	  vr = 0.0;
	  //vr = (-3.0*DISK_ALPHA*(1.0/Mach)*(1.0/Mach)*(1.0-delta_exp+xi_exp*pow(1./rs,-xi_exp)))*pow(r,2.);
        }else{
          omega = sqrt(O2);
          //vr = 0.0;
	  vr = -3.0/sqrt(r)*DISK_ALPHA*(1.0/Mach)*(1.0/Mach)*(1.0-delta_exp+xi_exp*pow(r/rs,-xi_exp));
        }





	
	//-------------------------
	//SET PRESSURE
	//------------------------
	if( rho<sim_RHO_FLOOR(theSim) ) rho = sim_RHO_FLOOR(theSim);
	double Pp = PoRho * rho;

	//---------------------------
	//R3B STUFF
	//---------------------------
	double xx = r*cos(phi);
	double yy = r*sin(phi);

	double u2 = sep*massratio/(1.+massratio);
	double u1 = sep - u2;


	double TU = xx*xx + yy*yy + 2.*(u1/dist_bh0 + u2/dist_bh1) ;
	double CJ = TU - pow((r*omega - r*pow(sep,-1.5)),2.); 
	

	double CJcrit = 100.0;
	double xL2 = 1.1;
	if (massratio == 1.0){
          CJcrit = 3.4568;
	  xL2 = 1.19841;
        }
        if (massratio == 0.1){
          CJcrit = 3.45154;
          xL2 = 1.25608;
        }
	if (massratio == 0.05){
	  CJcrit = 3.34669;
	  xL2 = 1.22557;
	}
        if (massratio == 0.03){
          CJcrit = 3.27402;
          xL2 =1.19963;
        }
	if (massratio == 0.01){
	  CJcrit = 3.15345;
	  xL2 = 1.14632;
	}
        if (massratio == 0.001){
          CJcrit = 3.03859;
	  xL2 = 1.06989;
        }


	double passive_scalar=0.000000001;
        //if (r< xL2 && CJ>CJcrit){
	if (r < xL2 && TU>CJcrit){
            passive_scalar = 1.;
            //printf("Setting Passive Scalar1");
        }

        double passive_scalar2=0.000000001;
        //if (r > xL2 && CJ>CJcrit){
	if (r > xL2 && TU>CJcrit){
            passive_scalar2 = 1.;
            //printf("Setting Passive Scalar1");
        }

	double passive_scalar3=0.000000001;
        if (CJ<CJcrit){
	  passive_scalar3 = 1.;
        }



	theCell->prim[RHO] = rho*exp(fac*(Pot0+Pot1)*Gam/cs2);
	theCell->prim[PPP] = Pp*exp(fac*(Pot0+Pot1)*Gam/cs2);
	theCell->prim[URR] = vr;
	theCell->prim[UPP] = omega - sim_rOm_a(theSim,r,sep)/r; //initial sep set to 1.
	theCell->prim[UZZ] = 0.0;
	theCell->prim[5] = passive_scalar;
	theCell->prim[6] = passive_scalar2;
	theCell->prim[7] = passive_scalar3;
	theCell->wiph = 0.0;
	theCell->divB = 0.0;
	theCell->GradPsi[0] = 0.0;
	theCell->GradPsi[1] = 0.0;
	theCell->GradPsi[2] = 0.0;
	//if(sim_NUM_C(theSim)<sim_NUM_Q(theSim)) theCell->prim[5] = passive_scalar;
	//if((sim_NUM_C(theSim)+1)<sim_NUM_Q(theSim)) theCell->prim[6] = passive_scalar2;
	//if((sim_NUM_C(theSim)+2)<sim_NUM_Q(theSim)) theCell->prim[7] = passive_scalar3;
}


void cell_init_MMRingCav(struct Cell ***theCells,struct Sim *theSim,struct MPIsetup * theMPIsetup) {
	
 // MM density params
  double rho0   = 1.0;
  double rs = 10.;
  double delta_exp   = 3.0;
  double xi_exp = 2.0;

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
  double phi0 = 0.0;
  double phi1 = M_PI;
  double eps0 = sim_G_EPS(theSim) * r0;
  double eps1 = sim_G_EPS(theSim) * r1;
  double eps = eps0/r0;

  if (massratio == 0.0){
    sep = 0.0;
    M1 = 0.0;
    r0 = 0.0;
    r1 = 0.0;
    eps0 = sim_G_EPS(theSim);
    eps1 = eps0;
  }


  double HoR = sim_HoR(theSim);   // for uniform Mach
  double Mach   = sqrt(Gam)*1./HoR;  
  double DISK_ALPHA = sim_EXPLICIT_VISCOSITY(theSim);

  int i, j,k;
  for (k = 0; k < sim_N(theSim,Z_DIR); k++) {
    for (i = 0; i < sim_N(theSim,R_DIR); i++) {
      double rm = sim_FacePos(theSim,i-1,R_DIR);
      double rp = sim_FacePos(theSim,i,R_DIR);
      double r = 0.5*(rm+rp);

      double xbh0 = r0;
      double xbh1 = -r1;
      double ybh0 = 0.0;
      double ybh1 = 0.0;

      double xx;
      double yy;
      double u1;
      double u2;
      double TU;
      double CJ;

      for (j = 0; j < sim_N_p(theSim,i); j++) {

        double phi = theCells[k][i][j].tiph - 0.5*theCells[k][i][j].dphi;
        double xpos = r*cos(phi);
        double ypos = r*sin(phi);

        double dist_bh0 = sqrt((xpos-xbh0)*(xpos-xbh0)+(ypos-ybh0)*(ypos-ybh0));
        double dist_bh1 = sqrt((xpos-xbh1)*(xpos-xbh1)+(ypos-ybh1)*(ypos-ybh1));


	//-------------------------
	//SET SOND SPEED FOR VERTICAL HYDROSTATIC EQLB for binary 
	//------------------------
	double OmegaK = 1./pow(r,1.5);
	double Omeff = sqrt(pow(dist_bh0*dist_bh0+eps*eps,-1.5)*M0+pow(dist_bh1*dist_bh1+eps*eps,-1.5)*M1);
	double Veff = sqrt(pow(dist_bh0*dist_bh0+eps*eps,-0.5)*M0+pow(dist_bh1*dist_bh1+eps*eps,-0.5)*M1); //Virial Thm

        //double cs = r*Omeff/Mach;
	double cs = Veff/Mach;
	double cs2 = cs*cs;
	double PoRho = cs2/Gam;


        double n = sim_PHI_ORDER(theSim);
        double Pot0 = M0/pow( pow(dist_bh0,n) + pow(eps0,n) , 1./n );
        double Pot1 = M1/pow( pow(dist_bh1,n) + pow(eps1,n) , 1./n );

	//-------------------------
	//FIND dnu/dr to set VR
	//------------------------
	double drcs2 = -1./(r*r)/Mach/Mach; // WRONG FOR Omeff
        double eps = sim_G_EPS(theSim);
        double disk_nu;
        double dr_nu;
        double DISK_ALPHA = sim_EXPLICIT_VISCOSITY(theSim);
        double dr_Omeff;
	double numer;
        if (sim_VISC_CONST(theSim)==1){
          disk_nu = DISK_ALPHA;
          dr_nu   = 0.0;
        }else{
          disk_nu = DISK_ALPHA*cs2/Omeff;
	  //
          numer = (  (r-r0 * cos(phi-phi0)) * pow(dist_bh0,-5.) )/(1.+massratio) + (  (r-r1 * cos(phi-phi1)) * pow(dist_bh1,-5.) )/(1.+1./massratio);
	  //
          dr_Omeff = - 1.5* numer / sqrt( pow(dist_bh0,-3.)/(1.+massratio) + pow(dist_bh1,-3.)/(1.+1./massratio)   );
	  //
	  dr_nu = DISK_ALPHA * ( drcs2/Omeff - cs2/Omeff/Omeff * dr_Omeff);
	}

	//-------------------------
	//SET DENS OMEGA
	//-------------------------	
	double rho = rho0*( pow(rs/r,delta_exp) )*exp(-pow(rs/r,xi_exp) );
        double omega = sqrt(1.)/pow(r,1.5);
        omega *= ( 1.+3./4.*(sep*sep/r/r*(massratio/((1.+massratio)*(1.+massratio)))) );
        double O2 = omega*omega + cs*cs/r/r*( 2.*rs*rs/r/r - 3. );
        double vr;
        if (r<3.0){
          omega = 1./pow(r,1.5);//sqrt(O2);
	  vr = 0.0;
          //vr = (-3.0*DISK_ALPHA*(1.0/Mach)*(1.0/Mach)*(1.0-delta_exp+xi_exp*pow(1./rs,-xi_exp)))*pow(r,2.);                                        
        }else{
          omega = sqrt(O2);
          //vr = 0.0;
	  vr = -3.0/sqrt(r)*DISK_ALPHA*(1.0/Mach)*(1.0/Mach)*(1.0-delta_exp+xi_exp*pow(r/rs,-xi_exp));
        }


	//-------------------------
	//SET PRESSURE
	//-------------------------
        if( rho<sim_RHO_FLOOR(theSim) ) rho = sim_RHO_FLOOR(theSim);
	double Pp = PoRho * rho;


	//--------------------------------
	// R3B stuff for Passive Scalar
	//-------------------------------
        xx = r*cos(phi);
        yy = r*sin(phi);

        u2 = sep*massratio/(1.+massratio);
        u1 = sep - u2;


        TU = xx*xx + yy*yy + 2.*(u1/dist_bh0 + u2/dist_bh1) ;
	CJ = TU - pow((r*omega - r*pow(sep,-1.5)),2);


        double CJcrit = 100.0;
	double xL2 = 1.1;
        if (massratio == 1.0){
          CJcrit = 3.4568;
          xL2 =1.19841;
	}
        if (massratio == 0.1){
          CJcrit = 3.45154;
          xL2 = 1.25608;
        }
        if (massratio == 0.05){
	  CJcrit = 3.34669;
          xL2 =1.22557;
	}
        if (massratio == 0.03){
          CJcrit = 3.27402;
          xL2 =1.19963;
        }
        if (massratio == 0.01){
          CJcrit = 3.15345;
          xL2 = 1.14632;
        }
        if (massratio == 0.001){
          CJcrit = 3.03859;
          xL2 =1.06989;
	}


	double passive_scalar=0.000000001;
        if (r< xL2 && TU>CJcrit){
	  passive_scalar = 1.;
	  //printf("Setting Passive Scalar1");                                                                                                     
        }

        double passive_scalar2=0.000000001;
        if (r > xL2 && TU>CJcrit){
	  passive_scalar2 = 1.;
	  //printf("Setting Passive Scalar1");                                                                                                     
        }

        double passive_scalar3=0.000000001;
        if (CJ<CJcrit){
          passive_scalar3 = 1.;
        }



        theCells[k][i][j].prim[RHO] = rho*exp(fac*(Pot0+Pot1)/cs2);
        theCells[k][i][j].prim[PPP] = Pp*exp(fac*(Pot0+Pot1)/cs2);
        theCells[k][i][j].prim[URR] = vr;
        theCells[k][i][j].prim[UPP] = omega -  sim_rOm_a(theSim,r,sep)/r; //initial sep set to 1. 
        theCells[k][i][j].prim[UZZ] = 0.0;
	theCells[k][i][j].prim[5] = passive_scalar;
	theCells[k][i][j].prim[6] = passive_scalar2;
	theCells[k][i][j].prim[7] = passive_scalar3;
        theCells[k][i][j].wiph = 0.0;
        theCells[k][i][j].divB = 0.0;
        theCells[k][i][j].GradPsi[0] = 0.0;
        theCells[k][i][j].GradPsi[1] = 0.0;
        theCells[k][i][j].GradPsi[2] = 0.0;
	//if((sim_NUM_C(theSim))<sim_NUM_Q(theSim)) theCells[k][i][j].prim[5] = passive_scalar;
	//if((sim_NUM_C(theSim)+1)<sim_NUM_Q(theSim)) theCells[k][i][j].prim[6] = passive_scalar2;
	//if((sim_NUM_C(theSim)+2)<sim_NUM_Q(theSim)) theCells[k][i][j].prim[7] = passive_scalar3;
      }
    }
  }

}
