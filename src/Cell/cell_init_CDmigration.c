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
	double phi0 = 0.0;
	double phi1 = M_PI; 
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

	
	double xbh0 = r0;
	double xbh1 = -r1;
	double ybh0 = 0.0;
	double ybh1 = 0.0;
			
				
	double phi = theCell->tiph-.5*theCell->dphi;
	double xpos = r*cos(phi);
	double ypos = r*sin(phi);
				
	double dist_bh0 = sqrt((xpos-xbh0)*(xpos-xbh0)+(ypos-ybh0)*(ypos-ybh0));
	double dist_bh1 = sqrt((xpos-xbh1)*(xpos-xbh1)+(ypos-ybh1)*(ypos-ybh1));
	

	double OmegaK = pow(r,-1.5);

	double cs = r*OmegaK/Mach;
	double cs2 = cs*cs;

	//cs0 = pow((dist_bh0*dist_bh0 + 0.05*0.05),-0.5) * sqrt(M0) * HoR;
	//cs1 = pow((dist_bh1*dist_bh1 + 0.05*0.05),-0.5) * sqrt(M1) * HoR;
	//PoRho0 =  cs0*cs0/Gam;
	//PoRho1 =  cs1*cs1/Gam; 
	
        //PoRho = PoRho0+PoRho1;
   	//PoRho = (cs0*cs0 + cs1*cs1)/Gam;
	double PoRho = cs2/Gam;   //for any gamma law eqn of state (adiabatic or iso)

				
	double n = sim_PHI_ORDER(theSim);
	double Pot0 = M0/pow( pow(dist_bh0,n) + pow(eps0,n) , 1./n );
	double Pot1 = M1/pow( pow(dist_bh1,n) + pow(eps1,n) , 1./n );
				


        //double cs    = sqrt(PoRho*Gam);//r*omega/Mach;
	//double scr_1 = sqrt(r1*r1 + r*r - 2.*r1*r*cos(phi1-phi));
	//double scr_0 = sqrt(r0*r0 + r*r - 2.*r0*r*cos(phi0-phi));
	//double drcs2 = - M0*(r-r0*cos(phi0 - phi))/pow((dist_bh0*dist_bh0 + 0.05*0.05), 1.5) - M1*(r-r1*cos(phi1 - phi))/pow((dist_bh1*dist_bh1 + 0.05*0.05), 1.5);


        //double cs2 = (cs0*cs0 + cs1*cs1);
	double drcs2 = -1./(r*r)/Mach/Mach;

 
	double eps = sim_G_EPS(theSim);
        double disk_nu;
	double dr_nu;
	double DISK_ALPHA = sim_EXPLICIT_VISCOSITY(theSim);
	double Omeff;
	double dr_Omeff;
	double numer;
        if (sim_VISC_CONST(theSim)==1){
          disk_nu = DISK_ALPHA;
	  dr_nu   = 0.0;
	}else{
	  Omeff = sqrt(pow(dist_bh0*dist_bh0+eps*eps,-1.5)*M0+pow(dist_bh1*dist_bh1+eps*eps,-1.5)*M1);
          disk_nu = DISK_ALPHA*cs2/Omeff;
	  //
	  numer = ( (r-r0 * cos(phi-phi0)) * pow(dist_bh0,-5.) )/(1.+massratio) + (  (r-r1 * cos(phi-phi1)) * pow(dist_bh1,-5.) )/(1.+1./massratio);
	  //
	  dr_Omeff = - 1.5* numer / sqrt( pow(dist_bh0,-3.)/(1.+massratio) + pow(dist_bh1,-3.)/(1.+1./massratio)   );
	  //
	  dr_nu = DISK_ALPHA * ( drcs2/Omeff - cs2/Omeff/Omeff * dr_Omeff);
	  //0.5*disk_nu/r; //0.5*DISK_ALPHA*pow(r,-0.5)/Mach/Mach;  //same as 0.5*disk_nu/r;
        }

	double omega;
	//if (MMOm==1){
	//omega = OmegaK*( 1.+3./4.*(sep*sep/r/r*(massratio/((1.+massratio)*(1.+massratio)))) );
	//double O2 = omega*omega + 1./2.*1/r * drcs2/Gam;    // 1/2 for rho -> sigma?- OmegaK*OmegaK/Mach/Mach/Gam;//add +1/(sigma r)*dP/dr trm 
	//}
	
	//if (TAddV==1){
	//double O2 = pow(dist_bh0*dist_bh0+eps*eps,-1.5)*M0 + pow(dist_bh1*dist_bh1+eps*eps,-1.5)*M1;
	//}
	//omega = sqrt(O2);

	//if (VirialV==1){
	double vp_vir = sqrt(pow(dist_bh0*dist_bh0+eps*eps,-0.5)*M0 + pow(dist_bh1*dist_bh1+eps*eps,-0.5)*M1);
	omega = vp_vir/r;
	//}



	double vr;
        double vrout;
        double rout = sim_RDAMP_OUTER(theSim);
	if (DISK_ALPHA > 0.0){
	  vr = -3./2.*disk_nu/r - 3.*dr_nu;
	  if ( sim_VISC_CONST(theSim)!=1 ){
              vrout = -3*DISK_ALPHA/Mach/Mach * pow(rout,-0.5);
              rho = rho/vr *vrout; // Choose rho so that we have constant Mdot in Damping region and rho =1.0 at rdamp_outer, less beyond
	  }
	}else{
	  vr = 0.0;
	}
	
	// st vr=0 inside a the cavity if it exists
	if (  r < (3. + 2.146*log10(massratio))   ) {//this function goes to 0 at q<0.04, =.9 at q=0.1, =3 at q=1
	  vr =0.0;
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
	double Pp = PoRho * rho;

	//R3B STUFF
	double xx = r*cos(phi);
	double yy = r*sin(phi);

	double u2 = sep*massratio/(1.+massratio);
	double u1 = sep - u2;


	double TU = xx*xx + yy*yy + 2.*(u1/dist_bh0 + u2/dist_bh1) ;
	double CJ = TU - pow((vp_vir - r*pow(sep,-1.5)),2.); 
	

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
	if (massratio == 0.01){
	  CJcrit = 3.15345;
	  xL2 = 1.14632;
	}
        if (massratio == 0.001){
          CJcrit = 3.03859;
	  xL2 = 1.06989;
        }


	double passive_scalar=0.000000001;
        if (r< xL2 && CJ>CJcrit){
            passive_scalar = 1.;
            //printf("Setting Passive Scalar1");
        }

        double passive_scalar2=0.000000001;
        if (r > xL2 && CJ>CJcrit){
            passive_scalar2 = 1.;
            //printf("Setting Passive Scalar1");
        }

	double passive_scalar3=0.000000001;
        if (CJ<CJcrit){
	  passive_scalar3 = 1.;
        }



	//        double passive_scalar2=0.000000001;
	// if (r>1.4 && r<1.5){
        //  if (phi<M_PI*(1./16.) && phi > -M_PI*(1./16.)){
        //    passive_scalar2 = 1.;
	//    //printf("Setting Passive Scalar1");
        //  }
        //}



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
  double phi0 = 0.0;
  double phi1 = M_PI;
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



	double OmegaK = 1./pow(r,1.5);

	double cs = r*OmegaK/Mach;
	double cs2 = cs*cs;
	double PoRho = cs2/Gam;


        double n = sim_PHI_ORDER(theSim);
        double Pot0 = M0/pow( pow(dist_bh0,n) + pow(eps0,n) , 1./n );
        double Pot1 = M1/pow( pow(dist_bh1,n) + pow(eps1,n) , 1./n );


	double rho = rho0;

       
	double drcs2 = -1./(r*r)/Mach/Mach; 


        double eps = sim_G_EPS(theSim);
        double disk_nu;
        double dr_nu;
        double DISK_ALPHA = sim_EXPLICIT_VISCOSITY(theSim);
        double Omeff;
        double dr_Omeff;
	double numer;
        if (sim_VISC_CONST(theSim)==1){
          disk_nu = DISK_ALPHA;
          dr_nu   = 0.0;
        }else{
          Omeff = sqrt(pow(dist_bh0*dist_bh0+eps*eps,-1.5)*M0+pow(dist_bh1*dist_bh1+eps*eps,-1.5)*M1);
	  //
          disk_nu = DISK_ALPHA*cs2/Omeff;
	  //
          numer = (  (r-r0 * cos(phi-phi0)) * pow(dist_bh0,-5.) )/(1.+massratio) + (  (r-r1 * cos(phi-phi1)) * pow(dist_bh1,-5.) )/(1.+1./massratio);
	  //
          dr_Omeff = - 1.5* numer / sqrt( pow(dist_bh0,-3.)/(1.+massratio) + pow(dist_bh1,-3.)/(1.+1./massratio)   );
	  //
	  dr_nu = DISK_ALPHA * ( drcs2/Omeff - cs2/Omeff/Omeff * dr_Omeff);
          //0.5*disk_nu/r; //0.5*DISK_ALPHA*pow(r,-0.5)/Mach/Mach;  //same as 0.5*disk_nu/r;
        }


	double omega;
        //if (MMOm==1){
        //omega = OmegaK*( 1.+3./4.*(sep*sep/r/r*(massratio/((1.+massratio)*(1.+massratio)))) );
        //double O2 = omega*omega + 1./2.*1/r * drcs2/Gam;    // 1/2 for rho -> sigma?- OmegaK*OmegaK/Mach/Mach/Gam; // add +1/(sigma r)*dP/dr tr 
	//	}

	// if (TAddV==1){
        //double O2 = pow(dist_bh0*dist_bh0+eps*eps,-1.5)*M0 + pow(dist_bh1*dist_bh1+eps*eps,-1.5)*M1;
        //}
	//omega = sqrt(O2);

	//if (VirialV==1){
	double vp_vir = sqrt(pow(dist_bh0*dist_bh0+eps*eps,-0.5)*M0 + pow(dist_bh1*dist_bh1+eps*eps,-0.5)*M1);
	omega = vp_vir/r; 
        //}  

	



        double vr;
	if (DISK_ALPHA > 0.0){
          vr = -3./2.*disk_nu/r - 3.*dr_nu;
        }else{
          vr = 0.0;
        }
	


      	  // st vr=0 inside a the cavity if it exists
	if (  r <  (3. + 2.146*log10(massratio))   ) {//this function goes to 0 at q<0.04, =1 at q=0.1, =3 at q=1
          vr =0.0;
	}





        if( rho<sim_RHO_FLOOR(theSim) ) rho = sim_RHO_FLOOR(theSim);
        
	 //double Pp = cs*cs * rho/Gam;
	double Pp = PoRho * rho;

        xx = r*cos(phi);
        yy = r*sin(phi);

        u2 = sep*massratio/(1.+massratio);
        u1 = sep - u2;


        TU = xx*xx + yy*yy + 2.*(u1/dist_bh0 + u2/dist_bh1) ;
	CJ = TU - pow((vp_vir - r*pow(sep,-1.5)),2);


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
        if (massratio == 0.01){
          CJcrit = 3.15345;
          xL2 = 1.14632;
        }
        if (massratio == 0.001){
          CJcrit = 3.03859;
          xL2 =1.06989;
	}


	double passive_scalar=0.000000001;
        if (r< xL2 && CJ>CJcrit){
	  passive_scalar = 1.;
	  //printf("Setting Passive Scalar1");                                                                                                     
        }

        double passive_scalar2=0.000000001;
        if (r > xL2 && CJ>CJcrit){
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
