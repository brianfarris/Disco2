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
  if (time_global <= 0.3 ){
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
    //Rot Frame or no here
    double tramp  = sim_tramp(theSim);
    if(mpisetup_MyProc(theMPIsetup)==0){
      char CstParamFilename[256];
      sprintf(CstParamFilename,"CstParam.dat");
      FILE * CstParamFile = fopen(CstParamFilename,"a");
      fprintf(CstParamFile,"%e %e %e %e %e %e %e %e %e %e %i %e %e %e \n", q, sep0, alpha, Mach0, MdoMs, rmin, rmax, tmigon, tacc, eps, nu_cst, Rcut, Gam, tramp);
      fclose(CstParamFile);
    }
  }

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

    double * TrVec_temp    = malloc(sizeof(double) * num_r_points_global); //DD: store Trqs each measure step 
    
    double TrScal_temp   = 0.0; // sum total trq each measure
    double TrScal_reduce = 0.0;

    double mass_near_bh0_r0p1_temp = 0.0;
    double mass_near_bh1_r0p1_temp = 0.0;
    double mass_near_bh0_r0p1_reduce = 0.0;
    double mass_near_bh1_r0p1_reduce = 0.0;

    double mass_near_bh0_r0p2_temp = 0.0;
    double mass_near_bh1_r0p2_temp = 0.0;
    double mass_near_bh0_r0p2_reduce = 0.0;
    double mass_near_bh1_r0p2_reduce = 0.0;

    double mass_near_bh0_r0p4_temp = 0.0;
    double mass_near_bh1_r0p4_temp = 0.0;
    double mass_near_bh0_r0p4_reduce = 0.0;
    double mass_near_bh1_r0p4_reduce = 0.0;


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
      TrVec_temp[i]=0.0;//DD - reset
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
          double phi = cell_tiph(cell_single(theCells,i,j,k));
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

	  double dist_bh0 = sqrt(r_bh0*r_bh0 + r*r - 2.*r_bh0*r*cos(phi_bh0-phi));
          double dist_bh1 = sqrt(r_bh1*r_bh1 + r*r - 2.*r_bh1*r*cos(phi_bh1-phi));

          // NOTE THIS IS THE derivative of the potential (add a -1/r for the force)                                                  
          double dPhi_dphi  =  abhbin*r/((1.+q)*(1.+1./q))*sin(phi-phi_bh0) * 
	    ( pow( (eps*eps + dist_bh0*dist_bh0),-1.5)- pow( (eps*eps + dist_bh1*dist_bh1) ,-1.5) ) ;
          //above DD: for general q - entire binary
          

          if (dist_bh0<0.1){
            double dM = rho*dV;
            mass_near_bh0_r0p1_temp +=dM;
          }
          if (dist_bh1<0.1){
            double dM = rho*dV;
            mass_near_bh1_r0p1_temp +=dM;
          }

          if (dist_bh0<0.2){
            double dM = rho*dV;
            mass_near_bh0_r0p2_temp +=dM;
          }
          if (dist_bh1<0.2){
            double dM = rho*dV;
            mass_near_bh1_r0p2_temp +=dM;
          }

          if (dist_bh0<0.4){
            double dM = rho*dV;
            mass_near_bh0_r0p4_temp +=dM;
          }
          if (dist_bh1<0.4){
            double dM = rho*dV;
            mass_near_bh1_r0p4_temp +=dM;
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
            EquatDiag_temp[(theDiagnostics->offset_eq+position)*NUM_EQ+15] = divB;
            EquatDiag_temp[(theDiagnostics->offset_eq+position)*NUM_EQ+16] = 180./M_PI*0.5*asin(-Br*Bp/(0.5*B2));
            EquatDiag_temp[(theDiagnostics->offset_eq+position)*NUM_EQ+17] = GradPsi_r;
	    EquatDiag_temp[(theDiagnostics->offset_eq+position)*NUM_EQ+18] = r*rho*dPhi_dphi/r; //DDwas 12
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
          VectorDiag_temp[(sim_N0(theSim,R_DIR)+i-imin)*NUM_VEC+8] += (vr*rho/sim_N_p(theSim,i)*dz);//(Bz/sim_N_p(theSim,i)*dz) ;
          VectorDiag_temp[(sim_N0(theSim,R_DIR)+i-imin)*NUM_VEC+9] += (E_hydro/sim_N_p(theSim,i)*dz) ;
          VectorDiag_temp[(sim_N0(theSim,R_DIR)+i-imin)*NUM_VEC+10] += (KE/sim_N_p(theSim,i)*dz) ;
          VectorDiag_temp[(sim_N0(theSim,R_DIR)+i-imin)*NUM_VEC+11] += (rhoe/sim_N_p(theSim,i)*dz) ;
          VectorDiag_temp[(sim_N0(theSim,R_DIR)+i-imin)*NUM_VEC+12] += (0.5*B2/sim_N_p(theSim,i)*dz) ;
          VectorDiag_temp[(sim_N0(theSim,R_DIR)+i-imin)*NUM_VEC+13] += (Br*Bp/sim_N_p(theSim,i)*dz) ;
          VectorDiag_temp[(sim_N0(theSim,R_DIR)+i-imin)*NUM_VEC+14] += (divB/sim_N_p(theSim,i)*dz) ;
          VectorDiag_temp[(sim_N0(theSim,R_DIR)+i-imin)*NUM_VEC+15] += (180./M_PI*0.5*asin(-Br*Bp/(0.5*B2))/sim_N_p(theSim,i)*dz) ;
          VectorDiag_temp[(sim_N0(theSim,R_DIR)+i-imin)*NUM_VEC+16] += (GradPsi_r/sim_N_p(theSim,i)*dz) ; 
	  // Don't count torques within a certain radius from secondary add r and 2pi for integration
	  if ((i>=imin_noghost) && (i<imax_noghost) && (k>=kmin_noghost) && (k<kmax_noghost)){
	    if (dist_bh1 > Rcut*Rhill && dist_bh0 > Rcut*Rhill){ //add 2pi to Trq to not azi average - see Scalar int below
	      VectorDiag_temp[(sim_N0(theSim,R_DIR)+i-imin)*NUM_VEC+17] += (2.*M_PI*r*rho*dPhi_dphi/sim_N_p(theSim,i)*dz); 
	    //positive gives torque ON binary                 
	      TrVec_temp[(sim_N0(theSim,R_DIR)+i-imin)]                 += (2.*M_PI*r*rho*dPhi_dphi/sim_N_p(theSim,i)*dz);//2pi/Nphi = dphi
	    }
	    }else{
	    VectorDiag_temp[(sim_N0(theSim,R_DIR)+i-imin)*NUM_VEC+17] += 0.0;
	    TrVec_temp[(sim_N0(theSim,R_DIR)+i-imin)]                 += 0.0;
          }
	  if (r<r_bh1){ ///Mass encolsed by binary               
            VectorDiag_temp[(sim_N0(theSim,R_DIR)+i-imin)*NUM_VEC+18] += (2*M_PI*rho/sim_N_p(theSim,i)*dz) ;
          }
          // the above are just placeholders. Put the real diagnostics you want here, then adjust NUM_DIAG accordingly.
        }
      }
    }


    for (i=imin;i<imax;++i){
      double rp = sim_FacePos(theSim,i,R_DIR);
      double rm = sim_FacePos(theSim,i-1,R_DIR);
      TrScal_temp += TrVec_temp[(sim_N0(theSim,R_DIR)+i-imin)]*(rp-rm);
      for (n=0;n<NUM_SCAL;++n){
        if (n==16){ // the vol. el. 'r' already added so dr not dr^2 - does not need extra 1/2 in /(dz) below
	  ScalarDiag_temp[n] += VectorDiag_temp[(sim_N0(theSim,R_DIR)+i-imin)*NUM_VEC+n+1]*(rp-rm);
        }else{
          // mult by delta r^2 because we are doing an r integration (needs extra factor of 1/2 comes with /(dz) below)
	  ScalarDiag_temp[n] += VectorDiag_temp[(sim_N0(theSim,R_DIR)+i-imin)*NUM_VEC+n+1]*(rp*rp-rm*rm);
        }

      }
    }



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
    MPI_Allreduce( &mass_near_bh1_r0p1_temp,&mass_near_bh1_r0p1_reduce , 1, MPI_DOUBLE, MPI_SUM, sim_comm);

    MPI_Allreduce( &mass_near_bh0_r0p2_temp,&mass_near_bh0_r0p2_reduce , 1, MPI_DOUBLE, MPI_SUM, sim_comm);
    MPI_Allreduce( &mass_near_bh1_r0p2_temp,&mass_near_bh1_r0p2_reduce , 1, MPI_DOUBLE, MPI_SUM, sim_comm);

    MPI_Allreduce( &mass_near_bh0_r0p4_temp,&mass_near_bh0_r0p4_reduce , 1, MPI_DOUBLE, MPI_SUM, sim_comm);
    MPI_Allreduce( &mass_near_bh1_r0p4_temp,&mass_near_bh1_r0p4_reduce , 1, MPI_DOUBLE, MPI_SUM, sim_comm);

    double mass_inside_1_reduce;
    MPI_Allreduce( &mass_inside_1_temp,&mass_inside_1_reduce , 1, MPI_DOUBLE, MPI_SUM, sim_comm);

    double mass_inner_wedge_reduce;
    MPI_Allreduce( &mass_inner_wedge_temp,&mass_inner_wedge_reduce , 1, MPI_DOUBLE, MPI_SUM, sim_comm);

    double RMIN = sim_MIN(theSim,R_DIR);
    double RMAX = sim_MAX(theSim,R_DIR);
    double ZMIN = sim_MIN(theSim,Z_DIR);
    double ZMAX = sim_MAX(theSim,Z_DIR);


    for (n=0;n<NUM_SCAL;++n){
      if (n!=16){ //Don't Vol avg Torques - only Z avg                     
        ScalarDiag_reduce[n] *= dtout/((ZMAX-ZMIN)*(RMAX*RMAX-RMIN*RMIN));
      }else if(n==16){
        ScalarDiag_reduce[n] *= dtout/(ZMAX-ZMIN);
	//TrScal_reduce        *= 1./(ZMAX-ZMIN);   don't volume average sigma = int{rho * dz}             
      }
    }


    int req1_found = 0;
    double Mdot_near_req1,r_near_req1;
    for (i=0;i<num_r_points_global;++i){
      double r = VectorDiag_reduce[i*NUM_VEC]/(ZMAX-ZMIN);
      if (r>1.0 && req1_found==0){
        Mdot_near_req1 = VectorDiag_reduce[i*NUM_VEC+5]*2.*M_PI*r/(ZMAX-ZMIN);
        r_near_req1 = r;
        //printf("req1_found i: %d, r: %e\n",i,r);
        req1_found = 1;
      }
      for (n=0;n<NUM_VEC;++n){
        VectorDiag_reduce[i*NUM_VEC+n] *= dtout/(ZMAX-ZMIN);
      }
    }

    if(mpisetup_MyProc(theMPIsetup)==0){
      char DiagMdotFilename[256];
      sprintf(DiagMdotFilename,"DiagMdot.dat");
      FILE * DiagMdotFile = fopen(DiagMdotFilename,"a");
      fprintf(DiagMdotFile,"%e %e %e %e %e %e %e %e %e %e %e %e %e\n",timestep_get_t(theTimeStep), Mdot_near_req1,r_near_req1,mass_inside_1_reduce,mass_inner_wedge_reduce,gravMass_Mdot(theGravMasses,0),gravMass_Mdot(theGravMasses,1),gravMass_Macc(theGravMasses,0),gravMass_Macc(theGravMasses,1),mass_near_bh0_r0p2_reduce,mass_near_bh1_r0p2_reduce,mass_near_bh0_r0p4_reduce,mass_near_bh1_r0p4_reduce);       
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
    double r_bh1 = gravMass_r(theGravMasses,1);
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

    double Mdotp = gravMass_Mdot(theGravMasses,0);
    double Maccp = gravMass_Macc(theGravMasses,0);
    double Mdots = gravMass_Mdot(theGravMasses,1);
    double Maccs = gravMass_Macc(theGravMasses,1);

    double sim_trq = gravMass_total_torque(theGravMasses,0);
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
      fprintf(DiagBpFile,"%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e \n",t, r_bh0, r_bh1, a_bin, phi_bh0, phi_bh1, ecc, E, Ltot, Mdotp, Maccp, Om, Om, sim_trq, TrScal_reduce, Mdots, Maccs);
      fclose(DiagBpFile);
    }
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
 
    free(TrVec_temp); //DD
    free(ScalarDiag_temp);
    free(ScalarDiag_reduce);
    free(VectorDiag_temp);
    free(VectorDiag_reduce);
    free(EquatDiag_temp);
    free(EquatDiag_reduce);

    theDiagnostics->tdiag_measure += theDiagnostics->dtdiag_measure;
  } 
}






