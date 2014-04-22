#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Sim.h"
#include "../Headers/Face.h"
#include "../Headers/GravMass.h"
#include "../Headers/header.h"

void cell_single_init_CDmigration(struct Cell *theCell, struct Sim *theSim,int i,int j,int k){

	double rho   = 1.0;
	
	
	double Gam = sim_GAMMALAW(theSim);
	double fac = 0.0001;

	double Mtotal = 1.0;
	double sep = sim_sep0(theSim);
	double massratio = sim_MassRatio(theSim);
	double M0 = Mtotal/(1.+massratio);
	double M1 = Mtotal/(1.+1./massratio);
	double r0 = sep*M1/Mtotal;
	double r1 = sep*M0/Mtotal;
	double eps0 = sim_G_EPS(theSim) * r0;
        double eps1 = sim_G_EPS(theSim) * r1;	


	if (massratio == 0.0){
	  sep = 0.0;
	  M1 = 0.0;
	  r0 = 0.0;
	  r1 = 0.0;	 
	  eps0 = sim_G_EPS(theSim);
	  eps1 = eps0;
	}

	//double Mach   = 10.0 *pow(sep,(-0.5)) ; //For uniform cs
	double HoR = sim_HoR(theSim);
	double Mach   = 1./HoR;                   //For uniform Mach	


	double rm = sim_FacePos(theSim,i-1,R_DIR);
	double rp = sim_FacePos(theSim,i,R_DIR);
	double r = 0.5*(rm+rp);

	// see below		
	//double cs;
	//cs = pow(r,(-0.5))/Mach;
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
				
	double dist_bh0 = sqrt((xpos-xbh0)*(xpos-xbh0)+(ypos-ybh0)*(ypos-ybh0));
	double dist_bh1 = sqrt((xpos-xbh1)*(xpos-xbh1)+(ypos-ybh1)*(ypos-ybh1));

	///```````````````````````````0
	double cs0,cs1;
        double PoRho0,PoRho1,PoRho;
        double vp0,vp1;
        if (dist_bh0>0.05){
	  vp0 = pow(dist_bh0,-0.5) * sqrt(M0) ;
	  cs0 = HoR*vp0;
	  PoRho0 = cs0*cs0/sim_GAMMALAW(theSim);
        } else{
	  vp0 = pow(0.05,-0.5);
	  cs0 = HoR*vp0;
	  PoRho0 = cs0*cs0/sim_GAMMALAW(theSim);
	  //PoRho = T0 * exp(-r*r/(r_a*r_a));                                                                                             
        }
        if (dist_bh1>0.05){
	  vp1 = pow(dist_bh1,-0.5) * sqrt(M1) ;
	  cs1 = HoR*vp1;
	  PoRho1 = cs1*cs1/sim_GAMMALAW(theSim);
        } else{
	  vp1 = pow(0.05,-0.5);
	  cs1 = HoR*vp1;
	  PoRho1 = cs1*cs1/sim_GAMMALAW(theSim);
	  //PoRho = T0 * exp(-r*r/(r_a*r_a));                                                                                             
        }

        PoRho = PoRho0+PoRho1;
	///

				
	double n = sim_PHI_ORDER(theSim);
	double Pot0 = M0/pow( pow(dist_bh0,n) + pow(eps0,n) , 1./n );
	double Pot1 = M1/pow( pow(dist_bh1,n) + pow(eps1,n) , 1./n );
				

	double omega = 1./pow(r,1.5);

        double cs = r*omega/Mach;
        //      double cs = pow(r,-0.5)/Mach;

 
	double eps = sim_G_EPS(theSim);
        double disk_nu;
	double dr_nu;
	double DISK_ALPHA = sim_EXPLICIT_VISCOSITY(theSim);
        if (sim_VISC_CONST(theSim)==1){
          disk_nu = DISK_ALPHA;
	  dr_nu   = 0.0;
	}else{
          disk_nu = DISK_ALPHA*cs*cs/sqrt(pow(dist_bh0*dist_bh0+eps*eps,-1.5)*M0+pow(dist_bh1*dist_bh1+eps*eps,-1.5)*M1);
	  dr_nu   = 0.5*DISK_ALPHA/Mach/Mach * sqrt(M0+M1) * pow(r,-0.5);
        }


	if (massratio!=0.0){
	  omega *= 1.+3./4.*(sep*sep/r/r*(massratio/((1.+massratio)*(1.+massratio))));
	}
	double O2 = omega*omega - omega*omega/Mach/Mach/Gam; // add 1/(sigma r) * dP/dr term
	omega = sqrt(O2);
   

	double vr;
	if (DISK_ALPHA > 0.0){
	  vr = -3./2.*disk_nu/r - 3.*dr_nu;
	}else{
	  vr = 0.0;
	}

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
				
	//double Pp = 1./20.*1./20.*rho/Gam; //cs*cs*rho/Gam * r;	//mult by r to get cst Pressure
	//double Pp = cs*cs*rho/Gam;
	double Pp = PoRho * rho;


	theCell->prim[RHO] = rho*exp(fac*(Pot0+Pot1)/cs/cs);
	theCell->prim[PPP] = Pp*exp(fac*(Pot0+Pot1)/cs/cs);
	theCell->prim[URR] = vr;
	theCell->prim[UPP] = omega - sim_rOm_a(theSim,r,sep)/r; //initial sep set to 1.
	theCell->prim[UZZ] = 0.0;
	theCell->wiph = 0.0;
	theCell->divB = 0.0;
	theCell->GradPsi[0] = 0.0;
	theCell->GradPsi[1] = 0.0;
	theCell->GradPsi[2] = 0.0;	
}

void cell_init_CDmigration(struct Cell ***theCells,struct Sim *theSim,struct MPIsetup * theMPIsetup) {
	
  double rho0   = 1.0;

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
  double eps0 = sim_G_EPS(theSim) * r0;
  double eps1 = sim_G_EPS(theSim) * r1;

  if (massratio == 0.0){
    sep = 0.0;
    M1 = 0.0;
    r0 = 0.0;
    r1 = 0.0;
    eps0 = sim_G_EPS(theSim);
    eps1 = eps0;
  }


 // double Mach   =  10.0 *pow(sep,(-0.5)) ; //For uniform cs
  double HoR = sim_HoR(theSim);   // for uniform Mach
  double Mach   = 1./HoR;
  
  double DISK_ALPHA = sim_EXPLICIT_VISCOSITY(theSim);

  int i, j,k;
  for (k = 0; k < sim_N(theSim,Z_DIR); k++) {
    for (i = 0; i < sim_N(theSim,R_DIR); i++) {
      double rm = sim_FacePos(theSim,i-1,R_DIR);
      double rp = sim_FacePos(theSim,i,R_DIR);
      double r = 0.5*(rm+rp);


      // See below
      //double cs;
      //cs = pow(r,(-0.5))/Mach;

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

        double dist_bh0 = sqrt((xpos-xbh0)*(xpos-xbh0)+(ypos-ybh0)*(ypos-ybh0));
        double dist_bh1 = sqrt((xpos-xbh1)*(xpos-xbh1)+(ypos-ybh1)*(ypos-ybh1));


        ///```````````````````````````0                                                                                                         

        double cs0,cs1;
        double PoRho0,PoRho1,PoRho;
        double vp0,vp1;
        if (dist_bh0>0.05){
          vp0 = pow(dist_bh0,-0.5) * sqrt(M0) ;
          cs0 = HoR*vp0;
          PoRho0 = cs0*cs0/sim_GAMMALAW(theSim);
        } else{
          vp0 = pow(0.05,-0.5);
          cs0 = HoR*vp0;
          PoRho0 = cs0*cs0/sim_GAMMALAW(theSim);
          //PoRho = T0 * exp(-r*r/(r_a*r_a));                                                                                                   
        }
        if (dist_bh1>0.05){
          vp1 = pow(dist_bh1,-0.5) * sqrt(M1) ;
          cs1 = HoR*vp1;
          PoRho1 = cs1*cs1/sim_GAMMALAW(theSim);
        } else{
          vp1 = pow(0.05,-0.5);
          cs1 = HoR*vp1;
          PoRho1 = cs1*cs1/sim_GAMMALAW(theSim);
          //PoRho = T0 * exp(-r*r/(r_a*r_a));                                                                                                   
        }

	PoRho = PoRho0+PoRho1;
	///  



        double n = sim_PHI_ORDER(theSim);
        double Pot0 = M0/pow( pow(dist_bh0,n) + pow(eps0,n) , 1./n );
        double Pot1 = M1/pow( pow(dist_bh1,n) + pow(eps1,n) , 1./n );

	    //---------------------------Shakura Sunyaev---------------------------//
		  //(7.21337*10^23 (1 - 0.0129692 Sqrt[1/R])^(7/10))/R^(3/4)  From Mathematica ".nb" 
		  // For 10^7Msun q=0.01 such that at r_sim=10 Mdisk = Msec 
	//double f = pow(1. - 0.0129692*sqrt(1./r),1./4.);
	double rho = rho0; //pow(r,(-3./4.)) * pow(f,(14./5.));
        double omega = 1./pow(r,1.5);
        double cs = r*omega/Mach;
        //        double cs = pow(r,-0.5)/Mach; 
	

	double eps = sim_G_EPS(theSim);
        double disk_nu;
	double dr_nu;
        double DISK_ALPHA = sim_EXPLICIT_VISCOSITY(theSim);
        if (sim_VISC_CONST(theSim)==1){
          disk_nu = DISK_ALPHA;
	  dr_nu   = 0.0;
	}else{
          disk_nu = DISK_ALPHA*cs*cs/sqrt(pow(dist_bh0*dist_bh0+eps*eps,-1.5)*M0+pow(dist_bh1*dist_bh1+eps*eps,-1.5)*M1);
          dr_nu   = 0.5*DISK_ALPHA/Mach/Mach * sqrt(M0+M1) * pow(r,-0.5);
        }



	if (massratio!=0.0){
	  omega *= 1.+3./4.*(sep*sep/r/r*(massratio/((1.+massratio)*(1.+massratio))));
	}
	
	double O2 = omega*omega - omega*omega/Mach/Mach/Gam; // 1/(r Sig)* dP/dr term

	omega = sqrt(O2);

        double vr;
        if (DISK_ALPHA > 0.0){
          vr = -3./2.*disk_nu/r - 3.*dr_nu;
        }else{
          vr = 0.0;
        }

        
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
        
        //double Pp = 1/20.*1/20.*rho/Gam; //cs*cs*rho/Gam * r; //mult by r to get cst P
	//double Pp = cs*cs*rho/Gam; 
	double Pp = PoRho * rho;

        theCells[k][i][j].prim[RHO] = rho*exp(fac*(Pot0+Pot1)/cs/cs);
        theCells[k][i][j].prim[PPP] = Pp*exp(fac*(Pot0+Pot1)/cs/cs);
        theCells[k][i][j].prim[URR] = vr;
        theCells[k][i][j].prim[UPP] = omega -  sim_rOm_a(theSim,r,sep)/r; //initial sep set to 1. 
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