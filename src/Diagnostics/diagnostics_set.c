#define DIAGNOSTICS_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Sim.h"
#include "../Headers/Cell.h"
#include "../Headers/GravMass.h"
#include "../Headers/Diagnostics.h"
#include "../Headers/TimeStep.h"
#include "../Headers/MPIsetup.h"
#include "../Headers/header.h"

void diagnostics_set(struct Diagnostics * theDiagnostics,struct Cell *** theCells,struct Sim * theSim,struct TimeStep * theTimeStep,struct MPIsetup * theMPIsetup,struct GravMass * theGravMasses){
  //  if ( (time_global <= 0.3 || sim_Restart(theSim)==1) && ii<2){
  if (iGlob<2){ 
  //Make a File that only writes 2-afew times  to output constant parameters    
    double q      = sim_MassRatio(theSim);
    double sep0   = sim_sep0(theSim);
    double alpha  = sim_EXPLICIT_VISCOSITY(theSim);
    double Mach0  = 1./sim_HoR(theSim);
    double MdoMs  = sim_Mdsk_o_Ms(theSim);
    double rmin   = sim_MIN(theSim, R_DIR);
    double rmax   = sim_MAX(theSim, R_DIR);
    double tmigon = sim_tmig_on(theSim);
    double tacc   = sim_RhoSinkOn(theSim);
    double eps    = sim_G_EPS(theSim);
    int nu_cst    = sim_VISC_CONST(theSim);    
    double Rcut   = sim_Rcut(theSim);
    double Gam    = sim_GAMMALAW(theSim);
    double Rsink0 = sim_Rsink0(theSim);
    double Rsink1 = sim_Rsink1(theSim);
    //Rot Frame or no here
    double tramp  = sim_tramp(theSim);
    double vRate  = sim_vRate(theSim);
    if(mpisetup_MyProc(theMPIsetup)==0){
      char CstParamFilename[256];
      sprintf(CstParamFilename,"CstParam.dat");
      FILE * CstParamFile = fopen(CstParamFilename,"a");
      fprintf(CstParamFile,"%e %e %e %e %e %e %e %e %e %e %i %e %e %e %e %e %e \n", q, sep0, alpha, Mach0, MdoMs, rmin, rmax, tmigon, tacc, eps, nu_cst, Rcut, Gam, tramp, vRate, Rsink0, Rsink1);
      fclose(CstParamFile);
    }
  }
  iGlob += 1;

 if (timestep_get_t(theTimeStep)>diagnostics_tdiag_measure(theDiagnostics)){
    int num_r_points = sim_N(theSim,R_DIR)-sim_Nghost_min(theSim,R_DIR)-sim_Nghost_max(theSim,R_DIR);
    int num_r_points_global = sim_N_global(theSim,R_DIR);

    int NUM_SCAL = theDiagnostics->NUM_DIAG; //DD: you can change the number of diagnostics in create and destroy
    int NUM_VEC = theDiagnostics->NUM_DIAG+1;
    int NUM_EQ = theDiagnostics->NUM_DIAG+2;

    int i,j,k,n;

    double * EquatDiag_temp = malloc(sizeof(double) * theDiagnostics->N_eq_cells*NUM_EQ);
    double * EquatDiag_reduce = malloc(sizeof(double) * theDiagnostics->N_eq_cells*NUM_EQ);
    double * VectorDiag_temp = malloc(sizeof(double) * num_r_points_global*NUM_VEC);
    double * VectorDiag_reduce = malloc(sizeof(double) * num_r_points_global*NUM_VEC);
    double * ScalarDiag_temp = malloc(sizeof(double)*NUM_SCAL);
    double * ScalarDiag_reduce = malloc(sizeof(double)*NUM_SCAL);
    //
    //double * TrVec_temp    = malloc(sizeof(double) * num_r_points_global); //DD: store Trqs each measure step 
    //
    double TrScal_temp   = 0.0; // sum total trq each measure
    double TrScal_reduce = 0.0;

    double mass_near_bh0_r0p1_temp = 0.0;
    double mass_near_bh1_0p5RH_temp = 0.0;
    double mass_near_bh0_r0p1_reduce = 0.0;
    double mass_near_bh1_0p5RH_reduce = 0.0;

    double mass_near_bh0_r0p2_temp = 0.0;
    double mass_near_bh1_RH_temp = 0.0;
    double mass_near_bh0_r0p2_reduce = 0.0;
    double mass_near_bh1_RH_reduce = 0.0;

    double mass_near_bh0_r0p4_temp = 0.0;
    double mass_near_bh1_2RH_temp = 0.0;
    double mass_near_bh0_r0p4_reduce = 0.0;
    double mass_near_bh1_2RH_reduce = 0.0;


    double dtout = timestep_get_t(theTimeStep)-theDiagnostics->toutprev;

    int imin = sim_Nghost_min(theSim,R_DIR);
    int imax = sim_N(theSim,R_DIR)-sim_Nghost_max(theSim,R_DIR);
    int kmin = sim_Nghost_min(theSim,Z_DIR);
    int kmax = sim_N(theSim,Z_DIR)-sim_Nghost_max(theSim,Z_DIR);


    int imin_noghost;
    if (sim_Nghost_min(theSim,R_DIR) == 1){
      imin_noghost = 0;
    } else{
      imin_noghost = sim_Nghost_min(theSim,R_DIR);
    }
    int kmin_noghost = sim_Nghost_min(theSim,Z_DIR);
    int imax_noghost = sim_N(theSim,R_DIR) - sim_Nghost_max(theSim,R_DIR);
    int kmax_noghost = sim_N(theSim,Z_DIR) - sim_Nghost_max(theSim,Z_DIR);



    for (n=0;n<NUM_SCAL;++n){
      ScalarDiag_temp[n]=0.0;
      ScalarDiag_reduce[n]=0.0;
    }
    for (i=0;i<num_r_points_global;++i){
      //TrVec_temp[i]=0.0;//DD - reset
      for (n=0;n<NUM_VEC;++n){
        VectorDiag_temp[i*NUM_VEC+n]=0.0;
      }
    }

    int position=0;
    for (i=0;i<num_r_points_global;++i){
      for(j = 0; j < theDiagnostics->N_p_global[i]; j++){
        for (n=0;n<NUM_EQ;++n){
          EquatDiag_temp[position]=0.0;
          ++position;
        }
      }
    }

    /*   // Fictitious forces on the Binary
    // FOR NOW ONLY FOR nu=cst and adv_arb movement 
    double Mbh0 = gravMass_M(theGravMasses,0);
    double Mbh1 = gravMass_M(theGravMasses,1);

    double rbh0 = gravMass_r(theGravMasses,0);
    double rbh1 = gravMass_r(theGravMasses,1);

    double a = rbh0 + rbh1;

    double nu = sim_EXPLICIT_VISCOSITY(theSim);
    double vr_sec = -sim_vRate(theSim)* 1.5 * nu/rbh1;
    double vr_prm = vr_sec * Mbh1/Mbh0;

    double wp_a = sim_rOm_a(theSim,rbh0,a);
    double rdrOmp_a = sim_rdrOm_a(theSim,rbh0,a);
    double dtOmp_a = sim_dtOm_a(theSim,rbh0,a);

    double ws_a = sim_rOm_a(theSim,rbh1,a);
    double rdrOms_a = sim_rdrOm_a(theSim,rbh1,a);
    double dtOms_a = sim_dtOm_a(theSim,rbh1,a);

    double Fsec_coriolis_phi = -2.*ws_a*vr_sec/rbh1;
    double Fsec_euler_phi =  -(rbh1*dtOms_a + vr_sec*rdrOms_a);
    double Fprm_coriolis_phi = -2.*wp_a*vr_prm/rbh0;
    double Fprm_euler_phi =  -(rbh0*dtOmp_a + vr_prm*rdrOmp_a);

    TrScal_temp += rbh0*(Fprm_coriolis_phi + Fprm_euler_phi)*Mbh0;
    TrScal_temp += rbh1*(Fsec_coriolis_phi + Fsec_euler_phi)*Mbh1;
    */

    double Rcut = sim_Rcut(theSim); //DD - for cutting out unresolved trq close to binary components
    double mass_inside_1_temp=0.0;
    position=0;
    for (k=kmin;k<kmax;++k){
      double zp = sim_FacePos(theSim,k,Z_DIR);
      double zm = sim_FacePos(theSim,k-1,Z_DIR);
      double z = 0.5*(zm+zp);
      double dz = zp-zm;
      for (i=imin;i<imax;++i){
        double rp = sim_FacePos(theSim,i,R_DIR);
        double rm = sim_FacePos(theSim,i-1,R_DIR);
        double r = 0.5*(rm+rp);
        double r_bh0 = gravMass_r(theGravMasses,0);
        double phi_bh0 = gravMass_phi(theGravMasses,0);
        double r_bh1 = gravMass_r(theGravMasses,1);
        double phi_bh1 = gravMass_phi(theGravMasses,1);
        double a = r_bh0 + r_bh1;        
        for (j=0;j<sim_N_p(theSim,i);++j){
	  //struct Cell * c = &(theCells[k][i][j]);
	  //double phi = c->tiph-.5*c->dphi;
	  //double dphi = c->dphi;
          double phi = cell_tiph(cell_single(theCells,i,j,k)) - 0.5*cell_dphi(cell_single(theCells,i,j,k));
          double dphi = cell_dphi(cell_single(theCells,i,j,k));
          double rho = cell_prim(cell_single(theCells,i,j,k),RHO);
          double press = cell_prim(cell_single(theCells,i,j,k),PPP);
          double vr = cell_prim(cell_single(theCells,i,j,k),URR);
          double vp_minus_w = cell_prim(cell_single(theCells,i,j,k),UPP)*r;
          double vp = cell_prim(cell_single(theCells,i,j,k),UPP)*r+sim_rOm_a(theSim,r,a);
          double vz = cell_prim(cell_single(theCells,i,j,k),UZZ);
          double Br = cell_prim(cell_single(theCells,i,j,k),BRR);
          double Bp = cell_prim(cell_single(theCells,i,j,k),BPP);
          double Bz = cell_prim(cell_single(theCells,i,j,k),BZZ);
          double v2   = vr*vr + vp*vp + vz*vz;
          double B2   = Br*Br + Bp*Bp + Bz*Bz;
          double rhoe = press/(sim_GAMMALAW(theSim)-1.);
          double KE = 0.5*rho*v2;
          double psi = cell_prim(cell_single(theCells,i,j,k),PSI);
          double E_hydro = rhoe + KE;
          double Cool = cell_Cool(theCells,i,j,k);
          double dV = 0.5*(rp*rp-rm*rm)*dphi;
          double divB = fabs(cell_divB(theCells,i,j,k)/dV);
          double GradPsi_r = cell_GradPsi(theCells,i,j,k,0)/dV;

          if (rp<=1.0){
            mass_inside_1_temp += rho*dV;
          }
          double passive_scalar = cell_prim(cell_single(theCells,i,j,k),5);
          double Omega = gravMass_omega(theGravMasses,1); //DD: changed from double Omega = 1.
	  double q = sim_MassRatio(theSim);       // Added Mass ratio DD 
          double t = timestep_get_t(theTimeStep);
          double eps = sim_G_EPS(theSim);
          
	  
	  //double abhbin = sqrt(r_bh0*r_bh0 + r_bh1*r_bh1 - 2.*r_bh0*r_bh1*cos(phi_bh1-phi_bh0)); //Added abhbin DD
	  double abhbin = r_bh0 + r_bh1;
	  double Rhill = abhbin * pow((q/3.),(1./3.)); //DD; Secondary Hill sphere

	  double xbh0 = r_bh0*cos(phi_bh0);
	  double ybh0 = r_bh0*sin(phi_bh0);

	  double xbh1 = r_bh1*cos(phi_bh1);
	  double ybh1 = r_bh1*sin(phi_bh1);

	  double xpos = r*cos(phi);
	  double ypos = r*sin(phi);

	  double dist_bh0 = sqrt((xpos-xbh0)*(xpos-xbh0)+(ypos-ybh0)*(ypos-ybh0));
	  double dist_bh1 = sqrt((xpos-xbh1)*(xpos-xbh1)+(ypos-ybh1)*(ypos-ybh1));

	  //double dist_bh0 = sqrt(r_bh0*r_bh0 + r*r - 2.*r_bh0*r*cos(phi_bh0-phi));
          //double dist_bh1 = sqrt(r_bh1*r_bh1 + r*r - 2.*r_bh1*r*cos(phi_bh1-phi));
          
          // NOTE THIS IS THE derivative of the potential (add a -1/r for the force)                                                  
          double dPhi_dphi  =  abhbin*r/((1.+q)*(1.+1./q))*sin(phi-phi_bh0) * 
	    ( pow( (eps*eps + dist_bh0*dist_bh0),-1.5)- pow( (eps*eps + dist_bh1*dist_bh1) ,-1.5) ) ;
          //above DD: for general q - entire binary
          double dPhi_dphi_s  =  -abhbin*r/((1.+q)*(1.+1./q))*sin(phi-phi_bh0) *
            (pow( (eps*eps + dist_bh1*dist_bh1) ,-1.5) ) ;
	  // above for secondary only
          //Now get the angel between phi direction and the lever arm to teh secondary to get torque only on secondary
	  double Cos_alph = (2.*r*r - 2.*r*abhbin/(1.+q) * cos(phi_bh1-phi))/dist_bh1;
	  double Sin_alph = sin(acos(Cos_alph));
	  

	  //In terms of force (for when eps!=0)
	  double fgrav( double M , double r , double eps, double n ){
	    return( M*pow(r,n-1.)/pow( pow(r,n) + pow(eps,n) ,1.+1./n) );
	  }

	  void gravMassForce( struct GravMass * theGravMasses ,struct Sim * theSim, int p , double r , double phi , double * fr , double * fp ){
	    double G_EPS=sim_G_EPS(theSim);
	    double PHI_ORDER=sim_PHI_ORDER(theSim);

	    double rp = gravMass_r(theGravMasses,p);
	    double pp = gravMass_phi(theGravMasses,p);
	    double cosp = cos(phi);
	    double sinp = sin(phi);
	    double dx = r*cosp-rp*cos(pp);
	    double dy = r*sinp-rp*sin(pp);
	    double script_r = sqrt(dx*dx+dy*dy);

	    double cosa = dx/script_r;
	    double sina = dy/script_r;

	    double cosap = cosa*cosp+sina*sinp;
	    double sinap = sina*cosp-cosa*sinp;

	    double f1;
	    if (sim_InitialDataType(theSim)==FIELDLOOP){
	      f1 = -fgrav_neg_centrifugal( gravMass_M(theGravMasses,p) , script_r , G_EPS, PHI_ORDER );
	    } else{
	      f1 = -fgrav( gravMass_M(theGravMasses,p) , script_r , G_EPS, PHI_ORDER );
	    }
	    *fr = cosap*f1;
	    *fp = sinap*f1;
	  }

	  double fr;
	  double fp;
	  double fBin_phi = 0.0;
	  double Fr = 0.0;
	  int p;
	  for( p=0 ; p<sim_NumGravMass(theSim); ++p ){
	    gravMassForce( theGravMasses , theSim , p , r , phi , &fr , &fp );
	    Fr += fr;
	    fBin_phi += fp;
	  }
	  
	  
	  //	  double fBin_phi = f0*sinap0 + f1*sinap1; //f1*sinap1;
	  //double fr0;
	  //double fp0;
	  //double fr1;
	  //double fp1;
	  //gravMassForce( theGravMasses , theSim , 0 , r , phi , &fr0 , &fp0 );
	  //gravMassForce( theGravMasses , theSim , 1 , r , phi , &fr1 , &fp1 );
          if ((i>=imin_noghost) && (i<imax_noghost) && (k>=kmin_noghost) && (k<kmax_noghost)){
            if (dist_bh1 > Rcut*Rhill){ //add 2pi to Trq to not azi average - see Scalar int below                                                  
	      //TrScal_temp -= ( r_bh0*fp0 * rho*r*dphi*(rp-rm) );
	      //TrScal_temp -= ( r_bh1*fp1 * rho*r*dphi*(rp-rm) );	      
              TrScal_temp   -= (r*rho* (fBin_phi) * r*dphi*(rp-rm));
            }
	  }	  	    

    
          if (dist_bh0<0.1){
            double dM = rho*dV;
            mass_near_bh0_r0p1_temp +=dM;
          }
          if (dist_bh1<0.5*Rhill){
            double dM = rho*dV;
            mass_near_bh1_0p5RH_temp +=dM;
          }

          if (dist_bh0<0.2){
            double dM = rho*dV;
            mass_near_bh0_r0p2_temp +=dM;
          }
          if (dist_bh1<Rhill){
            double dM = rho*dV;
            mass_near_bh1_RH_temp +=dM;
          }

          if (dist_bh0<0.4){
            double dM = rho*dV;
            mass_near_bh0_r0p4_temp +=dM;
          }
          if (dist_bh1<2.*Rhill){
            double dM = rho*dV;
            mass_near_bh1_2RH_temp +=dM;
          }



          if ((fabs(zp)<0.0000001)||(fabs(z)<0.0000001)){
            EquatDiag_temp[(theDiagnostics->offset_eq+position)*NUM_EQ+0] = r;
            EquatDiag_temp[(theDiagnostics->offset_eq+position)*NUM_EQ+1] = phi;
            EquatDiag_temp[(theDiagnostics->offset_eq+position)*NUM_EQ+2] = rho;
            EquatDiag_temp[(theDiagnostics->offset_eq+position)*NUM_EQ+3] = press;
            EquatDiag_temp[(theDiagnostics->offset_eq+position)*NUM_EQ+4] = vr;
            EquatDiag_temp[(theDiagnostics->offset_eq+position)*NUM_EQ+5] = vp;
            EquatDiag_temp[(theDiagnostics->offset_eq+position)*NUM_EQ+6] = vz;            
            EquatDiag_temp[(theDiagnostics->offset_eq+position)*NUM_EQ+7] = vp_minus_w;            
            EquatDiag_temp[(theDiagnostics->offset_eq+position)*NUM_EQ+8] = Cool;
            EquatDiag_temp[(theDiagnostics->offset_eq+position)*NUM_EQ+9] = Bz;            
            EquatDiag_temp[(theDiagnostics->offset_eq+position)*NUM_EQ+10] = E_hydro;            
            EquatDiag_temp[(theDiagnostics->offset_eq+position)*NUM_EQ+11] = KE;            
            EquatDiag_temp[(theDiagnostics->offset_eq+position)*NUM_EQ+12] = rhoe;
            EquatDiag_temp[(theDiagnostics->offset_eq+position)*NUM_EQ+13] = 0.5*B2;
            EquatDiag_temp[(theDiagnostics->offset_eq+position)*NUM_EQ+14] = Br*Bp;
            EquatDiag_temp[(theDiagnostics->offset_eq+position)*NUM_EQ+15] = dist_bh1*rho*dPhi_dphi_s/r * Sin_alph; //divB;
            EquatDiag_temp[(theDiagnostics->offset_eq+position)*NUM_EQ+16] = 180./M_PI*0.5*asin(-Br*Bp/(0.5*B2));
            EquatDiag_temp[(theDiagnostics->offset_eq+position)*NUM_EQ+17] = GradPsi_r;
	    EquatDiag_temp[(theDiagnostics->offset_eq+position)*NUM_EQ+18] = -r*rho*(fBin_phi);//-(r_bh0*fp0+r_bh1*fp1)*rho;
	                                                                                      //dPhi_dphi/r; //DDwas 12
            ++position;
          }
          // divide by number of phi cells to get phi average, mult by dz because we are doing a z integration;
          VectorDiag_temp[(sim_N0(theSim,R_DIR)+i-imin)*NUM_VEC+0] += (r/sim_N_p(theSim,i)*dz) ;
          VectorDiag_temp[(sim_N0(theSim,R_DIR)+i-imin)*NUM_VEC+1] += (rho/sim_N_p(theSim,i)*dz) ;
          VectorDiag_temp[(sim_N0(theSim,R_DIR)+i-imin)*NUM_VEC+2] += (press/sim_N_p(theSim,i)*dz) ;
          VectorDiag_temp[(sim_N0(theSim,R_DIR)+i-imin)*NUM_VEC+3] += (vr/sim_N_p(theSim,i)*dz) ;
          VectorDiag_temp[(sim_N0(theSim,R_DIR)+i-imin)*NUM_VEC+4] += (vp/sim_N_p(theSim,i)*dz) ;
          VectorDiag_temp[(sim_N0(theSim,R_DIR)+i-imin)*NUM_VEC+5] += (vz/sim_N_p(theSim,i)*dz) ;
          VectorDiag_temp[(sim_N0(theSim,R_DIR)+i-imin)*NUM_VEC+6] += (vp_minus_w/sim_N_p(theSim,i)*dz) ;
          VectorDiag_temp[(sim_N0(theSim,R_DIR)+i-imin)*NUM_VEC+7] += (Cool/sim_N_p(theSim,i)*dz) ;
          VectorDiag_temp[(sim_N0(theSim,R_DIR)+i-imin)*NUM_VEC+8] += (r * vr*rho/sim_N_p(theSim,i));//(Bz/sim_N_p(theSim,i)*dz) ;
          VectorDiag_temp[(sim_N0(theSim,R_DIR)+i-imin)*NUM_VEC+9] += (E_hydro/sim_N_p(theSim,i)*dz) ;
          VectorDiag_temp[(sim_N0(theSim,R_DIR)+i-imin)*NUM_VEC+10] += (KE/sim_N_p(theSim,i)*dz) ;
          VectorDiag_temp[(sim_N0(theSim,R_DIR)+i-imin)*NUM_VEC+11] += (rhoe/sim_N_p(theSim,i)*dz) ;
          VectorDiag_temp[(sim_N0(theSim,R_DIR)+i-imin)*NUM_VEC+12] += (0.5*B2/sim_N_p(theSim,i)*dz) ;
          VectorDiag_temp[(sim_N0(theSim,R_DIR)+i-imin)*NUM_VEC+13] += (Br*Bp/sim_N_p(theSim,i)*dz) ;
          VectorDiag_temp[(sim_N0(theSim,R_DIR)+i-imin)*NUM_VEC+14] += (divB/sim_N_p(theSim,i)*dz) ;
          VectorDiag_temp[(sim_N0(theSim,R_DIR)+i-imin)*NUM_VEC+15] += 0.0;//(180./M_PI*0.5*asin(-Br*Bp/(0.5*B2))/sim_N_p(theSim,i)*dz) ;
          VectorDiag_temp[(sim_N0(theSim,R_DIR)+i-imin)*NUM_VEC+16] += (GradPsi_r/sim_N_p(theSim,i)*dz) ; 
	  // Don't count torques within a certain radius from secondary add r and 2pi for integration
	  if ((i>=imin_noghost) && (i<imax_noghost) && (k>=kmin_noghost) && (k<kmax_noghost)){
	    if (dist_bh1 > Rcut*Rhill){ //add 2pi to Trq to not azi average - see Scalar int below
	      //VectorDiag_temp[(sim_N0(theSim,R_DIR)+i-imin)*NUM_VEC+17] += (2.*M_PI*r*rho*dPhi_dphi/sim_N_p(theSim,i)*dz);
	      //VectorDiag_temp[(sim_N0(theSim,R_DIR)+i-imin)*NUM_VEC+17] -= (r*2.*M_PI*rho*r*fBin_phi/sim_N_p(theSim,i));
	      VectorDiag_temp[(sim_N0(theSim,R_DIR)+i-imin)*NUM_VEC+17] -= (r*rho*fBin_phi * r*dphi);//((r_bh0*fp0+r_bh1*fp1)*rho *r*dphi);
	    //positive gives torque ON binary                 
	      //TrVec_temp[(sim_N0(theSim,R_DIR)+i-imin)]                 += (2.*M_PI*r*rho*dPhi_dphi/sim_N_p(theSim,i)*dz);//2pi/Nphi = dphi
	      //TrVec_temp[(sim_N0(theSim,R_DIR)+i-imin)]                 -= (r*2.*M_PI*rho*r*fBin_phi/sim_N_p(theSim,i));
	      //TrVec_temp[(sim_N0(theSim,R_DIR)+i-imin)]                 -= (r*rho*fBin_phi * r*dphi);
	      //TrScal_temp                                               -= (r*rho*fBin_phi * r*dphi*(rp-rm));
	    }
	    //}else{
	    // VectorDiag_temp[(sim_N0(theSim,R_DIR)+i-imin)*NUM_VEC+17] += 0.0;
	      //TrVec_temp[(sim_N0(theSim,R_DIR)+i-imin)]                 += 0.0;
          }
	  if (r<r_bh1){ ///Mass encolsed by binary               
            VectorDiag_temp[(sim_N0(theSim,R_DIR)+i-imin)*NUM_VEC+18] += ( 2*M_PI*rho/sim_N_p(theSim,i) );
          }
          // the above are just placeholders. Put the real diagnostics you want here, then adjust NUM_DIAG accordingly.
        }
      }
    }

    
     for (i=imin;i<imax;++i){
      for (j=0;j<sim_N_p(theSim,i);++j){
	k = 0;
	double rho = cell_prim(cell_single(theCells,i,j,k),RHO);
	VectorDiag_temp[(sim_N0(theSim,R_DIR)+i-imin)*NUM_VEC+15] += VectorDiag_temp[(sim_N0(theSim,R_DIR)+i-imin)*NUM_VEC+17]/(2*M_PI*rho)/sim_N_p(theSim,i);
      }    
     }


    for (i=imin;i<imax;++i){
      double rp = sim_FacePos(theSim,i,R_DIR);
      double rm = sim_FacePos(theSim,i-1,R_DIR);
      //TrScal_temp += TrVec_temp[(sim_N0(theSim,R_DIR)+i-imin)]*(rp-rm);
      for (n=0;n<NUM_SCAL;++n){
        if (n==16){ // the vol. el. 'r' already added so dr not dr^2 - does not need extra 1/2 in /(dz) below
	  ScalarDiag_temp[n] += VectorDiag_temp[(sim_N0(theSim,R_DIR)+i-imin)*NUM_VEC+n+1]*(rp-rm);
        }else{
          // mult by delta r^2 because we are doing an r integration r(rp-rm) = (rp+rm)/2*(rp-rm) = 1/2(rp^2-rm^2)
	  ScalarDiag_temp[n] += VectorDiag_temp[(sim_N0(theSim,R_DIR)+i-imin)*NUM_VEC+n+1]*(rp*rp-rm*rm);
        }

      }
    }

    // zero out vector array used to sum total Trq each timestep
    //for (i=imin;i<imax;++i){
    //  TrVec_temp[(sim_N0(theSim,R_DIR)+i-imin)] = 0.0;
    //}



    double mass_inner_wedge_temp = 0.0;
    if (mpisetup_MyProc(theMPIsetup)==0){
      //printf("myproc: %d, rm(0): %e, rp(0): %e\n",mpisetup_MyProc(theMPIsetup),sim_FacePos(theSim,-1,R_DIR),sim_FacePos(theSim,0,R_DIR));
      double rp = sim_FacePos(theSim,0,R_DIR);
      double rm = sim_FacePos(theSim,-1,R_DIR);
      for (j=0;j<sim_N_p(theSim,i);++j){
        double dphi = cell_dphi(cell_single(theCells,0,j,0));
        double dV = 0.5*(rp*rp-rm*rm)*dphi;
        double rho = cell_prim(cell_single(theCells,0,j,0),RHO);
        mass_inner_wedge_temp += rho*dV;
      }
    }
    MPI_Allreduce( ScalarDiag_temp,ScalarDiag_reduce , NUM_SCAL, MPI_DOUBLE, MPI_SUM, sim_comm);
    MPI_Allreduce( VectorDiag_temp,VectorDiag_reduce , num_r_points_global*NUM_VEC, MPI_DOUBLE, MPI_SUM, sim_comm);
    MPI_Allreduce( EquatDiag_temp,EquatDiag_reduce ,theDiagnostics->N_eq_cells*NUM_EQ, MPI_DOUBLE, MPI_SUM, sim_comm);

    MPI_Allreduce( &TrScal_temp, &TrScal_reduce , 1, MPI_DOUBLE, MPI_SUM, sim_comm); //DD

    MPI_Allreduce( &mass_near_bh0_r0p1_temp,&mass_near_bh0_r0p1_reduce , 1, MPI_DOUBLE, MPI_SUM, sim_comm);
    MPI_Allreduce( &mass_near_bh1_0p5RH_temp, &mass_near_bh1_0p5RH_reduce , 1, MPI_DOUBLE, MPI_SUM, sim_comm);

    MPI_Allreduce( &mass_near_bh0_r0p2_temp,&mass_near_bh0_r0p2_reduce , 1, MPI_DOUBLE, MPI_SUM, sim_comm);
    MPI_Allreduce( &mass_near_bh1_RH_temp,&mass_near_bh1_RH_reduce , 1, MPI_DOUBLE, MPI_SUM, sim_comm);

    MPI_Allreduce( &mass_near_bh0_r0p4_temp,&mass_near_bh0_r0p4_reduce , 1, MPI_DOUBLE, MPI_SUM, sim_comm);
    MPI_Allreduce( &mass_near_bh1_2RH_temp,&mass_near_bh1_2RH_reduce , 1, MPI_DOUBLE, MPI_SUM, sim_comm);

    double mass_inside_1_reduce;
    MPI_Allreduce( &mass_inside_1_temp,&mass_inside_1_reduce , 1, MPI_DOUBLE, MPI_SUM, sim_comm);

    double mass_inner_wedge_reduce;
    MPI_Allreduce( &mass_inner_wedge_temp,&mass_inner_wedge_reduce , 1, MPI_DOUBLE, MPI_SUM, sim_comm);

    double RMIN = sim_MIN(theSim,R_DIR);
    double RMAX = sim_MAX(theSim,R_DIR);
    double ZMIN = sim_MIN(theSim,Z_DIR);
    double ZMAX = sim_MAX(theSim,Z_DIR);


    for (n=0;n<NUM_SCAL;++n){
      if (n!=16 && n!=15){ //Don't Vol avg Torques - only Z avg                     
        ScalarDiag_reduce[n] *= dtout/((ZMAX-ZMIN)*(RMAX*RMAX-RMIN*RMIN));
	}else if(n==16 || n==15){
           ScalarDiag_reduce[n] *= dtout; ///(ZMAX-ZMIN);
	//TrScal_reduce        *= 1./(ZMAX-ZMIN);   don't volume average sigma = int{rho * dz}             
      }
    }

    double r_bh1 = gravMass_r(theGravMasses,1);
    int req1_found = 0;
    double Mdot_near_req1,r_near_req1;
    int reqs_found = 0;
    double Mdot_near_reqs,r_near_reqs;
    int reqOut_found = 0;
    double Mdot_near_reqOut,r_near_reqOut;
    int reqIn_found = 0;
    double Mdot_near_reqIn,r_near_reqIn;
    for (i=0;i<num_r_points_global;++i){
      double rr = VectorDiag_reduce[i*NUM_VEC]/(ZMAX-ZMIN);
      if (rr>1.0 && req1_found==0){
        Mdot_near_req1 = VectorDiag_reduce[i*NUM_VEC+8]*2.*M_PI*rr/(ZMAX-ZMIN);
        r_near_req1 = rr;
        //printf("req1_found i: %d, r: %e\n",i,r);
        req1_found = 1;
      }

      if (rr>r_bh1 && reqs_found==0){
	Mdot_near_reqs = VectorDiag_reduce[i*NUM_VEC+8]*2.*M_PI*rr/(ZMAX-ZMIN);
	r_near_reqs = rr;
	//printf("req1_found i: %d, r: %e\n",i,r);                                                                                              
	reqs_found = 1;
      }
      
      if (rr>r_bh1+0.5 && reqOut_found==0){
        Mdot_near_reqOut = VectorDiag_reduce[i*NUM_VEC+8]*2.*M_PI*rr/(ZMAX-ZMIN);
        r_near_reqOut = rr;
        //printf("req1_found i: %d, r: %e\n",i,r);                                                                                                 
        reqOut_found = 1;
      }

      if ((rr>r_bh1-0.5 || rr==sim_FacePos(theSim,i-1,R_DIR)) && reqIn_found==0){
        Mdot_near_reqIn = VectorDiag_reduce[i*NUM_VEC+8]*2.*M_PI*rr/(ZMAX-ZMIN);
        r_near_reqIn = rr;
        //printf("req1_found i: %d, r: %e\n",i,r); 
        reqIn_found = 1;
      }



      for (n=0;n<NUM_VEC;++n){
        VectorDiag_reduce[i*NUM_VEC+n] *= dtout/(ZMAX-ZMIN);
      }
    }

    if(mpisetup_MyProc(theMPIsetup)==0){
      char DiagMdotFilename[256];
      sprintf(DiagMdotFilename,"DiagMdot.dat");
      FILE * DiagMdotFile = fopen(DiagMdotFilename,"a");
      fprintf(DiagMdotFile,"%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",timestep_get_t(theTimeStep), Mdot_near_req1,r_near_req1, Mdot_near_reqs,r_near_reqs, Mdot_near_reqOut,r_near_reqOut, Mdot_near_reqIn,r_near_reqIn, mass_inside_1_reduce,mass_inner_wedge_reduce,gravMass_Mdot(theGravMasses,0),gravMass_Mdot(theGravMasses,1),gravMass_Macc(theGravMasses,0),gravMass_Macc(theGravMasses,1),mass_near_bh0_r0p2_reduce,mass_near_bh1_RH_reduce,mass_near_bh0_r0p4_reduce,mass_near_bh1_2RH_reduce,mass_near_bh0_r0p4_reduce,mass_near_bh1_0p5RH_reduce);       
      fclose(DiagMdotFile);

      char DiagTorqueFilename[256];
      sprintf(DiagTorqueFilename,"DiagTorque.dat");
      FILE * DiagTorqueFile = fopen(DiagTorqueFilename,"a");
      fprintf(DiagTorqueFile,"%e %e\n",timestep_get_t(theTimeStep), gravMass_total_torque(theGravMasses,0));       
      fclose(DiagTorqueFile);

    }


    //Binary params -DD ------ ADDED below DD---------//                      
    double t = timestep_get_t(theTimeStep);
    double r_bh0 = gravMass_r(theGravMasses,0);
    //double r_bh1 = gravMass_r(theGravMasses,1); defined above
    double phi_bh0 = gravMass_phi(theGravMasses,0);
    double phi_bh1 = gravMass_phi(theGravMasses,1);

    double Ltot = gravMass_Ltot(theGravMasses,1); //
    // double L1 = gravMass_L(theGravMasses,1); 
    double M0 = gravMass_M(theGravMasses,0);
    double M1 = gravMass_M(theGravMasses,1);
    double E  = gravMass_E(theGravMasses,1);
    double Om = gravMass_omega(theGravMasses,1);

    //double vr0 = gravMass_vr(theGravMasses,0);
    //double vr1 = gravMass_vr(theGravMasses,1);

    //double Mdotp = gravMass_Mdot(theGravMasses,0);
    //double Maccp = gravMass_Macc(theGravMasses,0);
    //double Mdots = gravMass_Mdot(theGravMasses,1);
    //double Maccs = gravMass_Macc(theGravMasses,1);

    //double sim_trq = gravMass_total_torque(theGravMasses,0);
    //If CM of binary strays from r=0, then phi1-phi0 != pi                     
    //double a_bin = sqrt(r_bh0*r_bh0 + r_bh1*r_bh1 - 2.*r_bh0*r_bh1*cos(phi_bh1-phi_bh0));
    double a_bin = r_bh0+r_bh1;
    double mu = M0*M1/(M0+M1);
    double l2 = mu*mu*pow(a_bin,4)*Om*Om;
    double ecc = sqrt( 1. + 2.*l2*E/(mu*(M0*M1)*(M0*M1)) );
    if(mpisetup_MyProc(theMPIsetup)==0){
      char DiagBPFilename[256];
      sprintf(DiagBPFilename,"BinaryParams.dat");
      FILE * DiagBpFile = fopen(DiagBPFilename,"a");
      fprintf(DiagBpFile,"%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e \n",t, r_bh0, r_bh1, a_bin, phi_bh0, phi_bh1, ecc, E, Ltot, Om, Om, Om, Om, gravMass_total_torque(theGravMasses,0), TrScal_reduce);
      fclose(DiagBpFile);
    }

    TrScal_reduce = 0.0;
    TrScal_temp = 0.0;
    //----------------ADDED above DD-------------------------//                                                                                                                



    //We are doing time averaged diagnostics, so mult by delta t and add it
    //We will divide by the total delta next time we save to disk;
    for (i=0;i<num_r_points_global;++i){
      for (n=0;n<NUM_VEC;++n){
        theDiagnostics->VectorDiag[i][n] += VectorDiag_reduce[i*NUM_VEC+n] ;
      }
    }
    for (n=0;n<NUM_SCAL;++n){
      theDiagnostics->ScalarDiag[n] += ScalarDiag_reduce[n] ;
    }

    position=0;
    for (i=0;i<num_r_points_global;++i){
      for(j = 0; j < theDiagnostics->N_p_global[i]; j++){
        for (n=0;n<NUM_EQ;++n){
          theDiagnostics->EquatDiag[position][n] = EquatDiag_reduce[position*NUM_EQ+n];
        }
        ++position;
      }
    }

    //update output time;
    theDiagnostics->toutprev = timestep_get_t(theTimeStep);
 
    //free(TrVec_temp); //DD
    free(ScalarDiag_temp);
    free(ScalarDiag_reduce);
    free(VectorDiag_temp);
    free(VectorDiag_reduce);
    free(EquatDiag_temp);
    free(EquatDiag_reduce);

    theDiagnostics->tdiag_measure += theDiagnostics->dtdiag_measure;
  } 
}






