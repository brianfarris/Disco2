#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Sim.h"
#include "../Headers/Face.h"
#include "../Headers/GravMass.h"
#include "../Headers/header.h"

double fgrav( double M , double r , double eps, double n ){
  return( M*pow(r,n-1.)/pow( pow(r,n) + pow(eps,n) ,1.+1./n) );
}

double fgrav_neg_centrifugal( double M , double r , double eps, double n ){
  double Om = 0.0;
  return( M*r*Om*Om );
}


void get_rho_sink( struct GravMass * theGravMasses, struct Sim * theSim, int p, double dt, double r0,double r1,double M0, double M1, double rho, double P, double * drho_dt_sink){

  double Mtotal = M0+M1;
  double alpha = sim_EXPLICIT_VISCOSITY(theSim);
  if (alpha <= 0.0){
    alpha=0.000025;
  }

  double one_o_nu = 1.0/(alpha*P/rho)*sqrt(pow(r0,-3)*M0/Mtotal + pow(r1,-3)*M1/Mtotal);
  //double one_o_nu1 = 1.0/(alpha*P/rho)*sqrt(pow(r1,-3)*M1/Mtotal);

  if (sim_VISC_CONST(theSim)==1){
    //alpha = sim_EXPLICIT_VISCOSITY(theSim);// / sim_HoR(theSim)/sim_HoR(theSim);
    one_o_nu = 1./sim_EXPLICIT_VISCOSITY(theSim);
  }

  double t_visc0 = TVISC_FAC * 2./3. * r0*r0*one_o_nu;
  double t_visc1 = TVISC_FAC * 2./3. * r1*r1*one_o_nu;


  if (t_visc0 < 10.* dt){
    t_visc0 = 10.*dt;
  }
  if (t_visc1 < 10.* dt){
    t_visc1 = 10.*dt;
  }

  *drho_dt_sink = 0.0;
  double rbh0 = gravMass_r(theGravMasses,0);
  double rbh1 = gravMass_r(theGravMasses,1);

  double a = rbh0+rbh1;


  //double Rbondi = 2.*gravMass_M(theGravMasses,p)/ (P/rho);
  double Rhill = a * pow(sim_MassRatio(theSim)/3., (1./3.));
  double Rsink1 = fmin(Rhill, sim_Rsink1(theSim));

  // double sink_size = fmin(Rhill, Rbondi);
  //sink_size = 0.5;
  //double sink_size = 1.0*Rhill; //0.5;
  if (a<=0.0){
    //sink_size = 0.05;
    sim_set_Rsink0(theSim, 0.05);
    sim_set_Rsink1(theSim, 0.05);
  }


  if (rho > sim_RHO_FLOOR(theSim)){
    if (p==0){
      //      if (r0<sim_Rsink0(theSim)  ){
      if (r0<sim_Rsink0(theSim)){
	*drho_dt_sink = rho / t_visc0;
      }
    } else if (p==1){
      //if (r1<sim_Rsink1(theSim)){
      if (r1<Rsink1){
	*drho_dt_sink = rho / t_visc1;
      }

    } else{
      exit(1);
    }
  }
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
  //double rH = theGravMasses[p].r*pow(theGravMasses[p].M/theGravMasses[0].M/3.,1./3.);
  //double pd = 0.8;
  //double fd = 1./(1.+exp(-( script_r/rH-pd)/(pd/10.)));

  //*fr = cosap*f1*fd;
  //*fp = sinap*f1*fd;
  *fr = cosap*f1;
  *fp = sinap*f1;
}

void cell_add_src( struct Cell *** theCells ,struct Sim * theSim, struct GravMass * theGravMasses , double dt ){
  ///Set a density scale for feedback to holes                                                                                                
  double Rcut = sim_Rcut(theSim);
  double dens_scale = ( sim_Mdsk_o_Ms(theSim)/( 1.+1./sim_MassRatio(theSim) ) )/( M_PI*sim_sep0(theSim)*sim_sep0(theSim) );
  if (   time_global <=    2.*M_PI*(sim_tmig_on(theSim) + sim_tramp(theSim))   ){
    dens_scale *= (time_global - 2*M_PI*sim_tmig_on(theSim))/( sim_tramp(theSim)*2*M_PI );
    // ABOVE^ Allow migration to turn on slowly                                                                                               
  }
  
  if (dens_scale <= 0.0){
    dens_scale = 0.0;
  }
    
  //dens_scale=1.0; //TESTING TESTING ************

  int GRAV2D=sim_GRAV2D(theSim);
  int POWELL=sim_POWELL(theSim);
  int i,j,k;
  double Mdot_temp[2];
  Mdot_temp[0] = 0.0;
  Mdot_temp[1] = 0.0;

  int imin_noghost;
  if (sim_Nghost_min(theSim,R_DIR) == 1){
    imin_noghost = 0;
  } else{
    imin_noghost = sim_Nghost_min(theSim,R_DIR);
  }
  int kmin_noghost = sim_Nghost_min(theSim,Z_DIR);
  int imax_noghost = sim_N(theSim,R_DIR) - sim_Nghost_max(theSim,R_DIR);
  int kmax_noghost = sim_N(theSim,Z_DIR) - sim_Nghost_max(theSim,Z_DIR);

  double ttorque_temp = 0.0;
  //  double total_Fr_temp = 0.0;


  double Mtotal = 1.0;
  //double sep = 1.0;
  double M0 = gravMass_M(theGravMasses,0);
  double M1 = gravMass_M(theGravMasses,1);

  double rbh0 = gravMass_r(theGravMasses,0);
  double rbh1 = gravMass_r(theGravMasses,1);
  
  double a = rbh0+rbh1;

  double phi_bh0 = gravMass_phi(theGravMasses,0);
  double phi_bh1 = gravMass_phi(theGravMasses,1);

  double xbh0 = rbh0*cos(phi_bh0);
  double ybh0 = rbh0*sin(phi_bh0);

  double xbh1 = rbh1*cos(phi_bh1);
  double ybh1 = rbh1*sin(phi_bh1);

  /*
  // FOR NOW ONLY FOR nu=cst and adv_arb movement                                                                                              
  double nu = sim_EXPLICIT_VISCOSITY(theSim);
  double vr_sec = -sim_vRate(theSim) * 1.5 * nu/rbh1;
  double vr_prm = vr_sec * M1/M0;

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

  ttorque_temp += rbh0*(Fprm_coriolis_phi + Fprm_euler_phi)*M0;
  ttorque_temp += rbh1*(Fsec_coriolis_phi + Fsec_euler_phi)*M1;
  */


  for( k=0 ; k<sim_N(theSim,Z_DIR) ; ++k ){
    double zm = sim_FacePos(theSim,k-1,Z_DIR);
    double zp = sim_FacePos(theSim,k,Z_DIR);
    for( i=0 ; i<sim_N(theSim,R_DIR) ; ++i ){
      double rm = sim_FacePos(theSim,i-1,R_DIR);
      double rp = sim_FacePos(theSim,i,R_DIR);
      for( j=0 ; j<sim_N_p(theSim,i) ; ++j ){
        struct Cell * c = &(theCells[k][i][j]);
        double phi = c->tiph-.5*c->dphi;
        double dphi = c->dphi;

        double rho = c->prim[RHO];
        double Pp  = c->prim[PPP];
        double r   = .5*(rp+rm);
        double vr  = c->prim[URR];
        double vz  = c->prim[UZZ];
        double vp  = c->prim[UPP]*r;
        double dz = zp-zm;
        //double dV = dphi*.5*(rp*rp-rm*rm)*dz;
	double dV = (rp-rm)*r*dphi*dz; //-DD

        double xpos = r*cos(phi);
        double ypos = r*sin(phi);

        double r0 = sqrt((xpos-xbh0)*(xpos-xbh0)+(ypos-ybh0)*(ypos-ybh0));
        double r1 = sqrt((xpos-xbh1)*(xpos-xbh1)+(ypos-ybh1)*(ypos-ybh1));


	double Rhill = a * pow(sim_MassRatio(theSim)/3., (1./3.));

        double z  = .5*(zp+zm);
        double R  = sqrt(r*r+z*z);
        double sint;
        double cost;
        double gravdist;
        if (GRAV2D==1){
          sint = 1.0;
          cost = 0.0;
          gravdist = r;
        }else{
          sint = r/R;
          cost = z/R;
          gravdist = R;
        }
        double fr,fp;
        double Fr = 0.0;
        double Fp = 0.0;
	//double Fmr =0.0; // for torques
	//double Fmp = 0.0;

	double dm = dens_scale * fabs(rho*dV/dz); // for torques - note: dV contains a factor of dz=2
 	
        int p;
        for( p=0 ; p<sim_NumGravMass(theSim); ++p ){
          gravMassForce( theGravMasses , theSim , p , gravdist , phi , &fr , &fp );
          Fr += fr;
          Fp += fp;
	  //Fmr -= fr*dm;
	  //Fmp = -fp*dm;
        }
        double w_a = sim_rOm_a(theSim,r,a);
        double rdrOm_a = sim_rdrOm_a(theSim,r,a);
        double F_centrifugal_r = w_a*w_a/r;
        double F_coriolis_r =  2.*w_a*vp/r;
        double F_coriolis_phi = -2.*w_a*vr/r;
        double dtOm_a = sim_dtOm_a(theSim,r,a);
        double F_euler_phi =  -(r*dtOm_a + vr*rdrOm_a);


        c->cons[SRR] += dt*dV*( rho*vp*vp + Pp )/r;
        c->cons[SRR] += dt*dV*rho*(Fr*sint+F_centrifugal_r/*+F_coriolis_r*/); //DD these terms are in split solver
	//c->cons[SRR] += dt*dV*rho*(Fr*sint+F_centrifugal_r +F_coriolis_r);
        c->cons[LLL] += dt*dV*rho*(Fp /*+F_coriolis_phi+F_euler_phi*/ -r*dtOm_a)*r; // these terms are in split
	//c->cons[LLL] += dt*dV*rho*(Fp + F_coriolis_phi+F_euler_phi /*-r*dtOm_a*/)*r;
        c->cons[SZZ] += dt*dV*rho*Fr*cost;
        c->cons[TAU] += dt*dV*rho*( (Fr*sint+F_centrifugal_r)*vr+Fr*vz*cost + (Fp+F_euler_phi)*vp);
	//double fpsec;
	//double frsec;
	//double fr0;
	//double fp0;
	//double fr1;
	//double fp1;
        if ((i>=imin_noghost) && (i<imax_noghost) && (k>=kmin_noghost) && (k<kmax_noghost)){
	  if (r1 > Rcut*Rhill){
	    //gravMassForce( theGravMasses , theSim , 0 , gravdist , phi , &fr0 , &fp0 );
	    //gravMassForce( theGravMasses , theSim , 1 , gravdist , phi , &fr1 , &fp1 );
	    ttorque_temp -= r*Fp*dm;
	    //ttorque_temp -= rbh0*fp0*dm;
	    //ttorque_temp -= rbh1*fp1*dm;
	  }
        }

	if (sim_RhoSinkOn(theSim)==1){
          double drho_dt_sink;
	  int num_sinks;
	  if (sim_GravMassType(theSim) <= 2){
	    num_sinks = sim_GravMassType(theSim);
	  }else{
	    num_sinks = 2;
	  }
          for( p=0 ; p<num_sinks; ++p ){
            get_rho_sink( theGravMasses, theSim, p, dt, r0, r1, M0, M1, rho, Pp, &drho_dt_sink);
            c->cons[RHO] -= drho_dt_sink * dt * dV;
            if ((i>=imin_noghost) && (i<imax_noghost) && (k>=kmin_noghost) && (k<kmax_noghost)){
              Mdot_temp[p] += drho_dt_sink*dV;
            }
          }
	}

        if(sim_runtype(theSim)==MHD){
          double Br = c->prim[BRR];
          double Bp = c->prim[BPP];
          double Bz = c->prim[BZZ];
          double B2 = Br*Br+Bp*Bp+Bz*Bz;
          double divB = c->divB;
          double GradPsi[3];
          GradPsi[0] = r*c->GradPsi[0];
          GradPsi[1] = r*c->GradPsi[1];
          GradPsi[2] = r*c->GradPsi[2];

          double vdotB = vr*Br+(vp/*+pow(r,-.5)*/)*Bp+vz*Bz;
          double BdotGradPsi = /*dV**/(Br*GradPsi[0]+Bp*GradPsi[1]+Bz*GradPsi[2]);
          double vdotGradPsi = /*dV**/(vr*GradPsi[0]+(vp-w_a)*GradPsi[1]+vz*GradPsi[2]);

          c->cons[SRR] += dt*dV*( .5*B2 - Bp*Bp )/r;
          c->cons[TAU] += dt*dV*Br*Bp*rdrOm_a;
          double psi = c->prim[PSI];
          c->cons[BRR] += dt*dV*( psi )/r;
          //cons[PSI] -= dt*dV*psi*DIVB_CH/DIVB_L;//DIVB_CH2/DIVB_CP2;
          //c->cons[BPP] += dt*dV*Br*sim_OM_A_DERIV(theSim,r);
          c->cons[BPP] += dt*dV*Br*rdrOm_a/r;

          c->cons[SRR] -= POWELL*dt*divB*Br;
          c->cons[SZZ] -= POWELL*dt*divB*Bz;
          c->cons[LLL] -= POWELL*dt*divB*Bp*r;
          c->cons[TAU] -= POWELL*dt*(divB*vdotB+BdotGradPsi);
          c->cons[BRR] -= POWELL*dt*divB*vr;
          c->cons[BPP] -= POWELL*dt*divB*vp/r;
          c->cons[BZZ] -= POWELL*dt*divB*vz;
          c->cons[PSI] -= POWELL*dt*vdotGradPsi;

        } 
        //cooling
        double numfac = 1.;
        double a_o_M = 100.0;
        if (sim_COOLING(theSim)==1){
          double alpha = sim_EXPLICIT_VISCOSITY(theSim);
          double cooling_cap = c->cons[TAU]*1000000000000000.0;
          //double cooling_term = numfac * pow(Pp/rho,0.5)*pow(r,-1.5)*pow(a_o_M,0.5)*dt*dV;
          double cooling_term = 9./4.*alpha*pow(PoRho_r1,-3)/rho*pow(Pp/rho,4.)*dt*dV;
          if (cooling_term>cooling_cap) cooling_term = cooling_cap;
          c->cons[TAU] -= cooling_term;
          c->Cool = cooling_term/dt/dV;
        }
      }
    }
  }



  double Mdot_reduce[2];
  MPI_Allreduce( Mdot_temp, Mdot_reduce , 2, MPI_DOUBLE, MPI_SUM, sim_comm);
  double ttorque_reduce;
  MPI_Allreduce(&ttorque_temp, &ttorque_reduce, 1, MPI_DOUBLE, MPI_SUM, sim_comm);
  int p;
  for( p=0 ; p<sim_NumGravMass(theSim); ++p ){
    gravMass_set_Mdot(theGravMasses,Mdot_reduce[p],p);  //dont set for now DD
    double Macc_prev = gravMass_Macc(theGravMasses,p);
    gravMass_set_Macc(theGravMasses,Macc_prev+Mdot_reduce[p]*dt,p);
    gravMass_set_total_torque(theGravMasses,ttorque_reduce,p);
  }
}

void cell_add_visc_src( struct Cell *** theCells ,struct Sim * theSim, struct GravMass * theGravMasses, double dt ){
  double Mtotal = 1.0;
  //double sep = 1.0;
  double M0 = gravMass_M(theGravMasses,0);
  double M1 = gravMass_M(theGravMasses,1);

  double r0 = gravMass_r(theGravMasses,0);
  double r1 = gravMass_r(theGravMasses,1);

  double a = r0 + r1;

  double phi_bh0 = gravMass_phi(theGravMasses,0);
  double phi_bh1 = gravMass_phi(theGravMasses,1);

  double xbh0 = r0*cos(phi_bh0);
  double ybh0 = r0*sin(phi_bh0);

  double xbh1 = r1*cos(phi_bh1);
  double ybh1 = r1*sin(phi_bh1);


  int i,j,k;
  for( k=0 ; k<sim_N(theSim,Z_DIR) ; ++k ){
    double zm = sim_FacePos(theSim,k-1,Z_DIR);
    double zp = sim_FacePos(theSim,k,Z_DIR);
    for( i=0 ; i<sim_N(theSim,R_DIR) ; ++i ){
      double rm = sim_FacePos(theSim,i-1,R_DIR);
      double rp = sim_FacePos(theSim,i,R_DIR);
      for( j=0 ; j<sim_N_p(theSim,i) ; ++j ){
        struct Cell * c = &(theCells[k][i][j]);
        double dphi = c->dphi;
        double rho = c->prim[RHO];
        double P = c->prim[PPP];
        double r   = .5*(rp+rm);
        double vr  = c->prim[URR];
        double dz = zp-zm;
        double dV = dphi*.5*(rp*rp-rm*rm)*dz;

        double phi = c->tiph - 0.5*c->dphi;

        double xpos = r*cos(phi);
        double ypos = r*sin(phi);

        double dist_bh0 = sqrt((xpos-xbh0)*(xpos-xbh0)+(ypos-ybh0)*(ypos-ybh0));
        double dist_bh1 = sqrt((xpos-xbh1)*(xpos-xbh1)+(ypos-ybh1)*(ypos-ybh1));



        double alpha = sim_EXPLICIT_VISCOSITY(theSim);
        double rdrOm_a = sim_rdrOm_a(theSim,r,a);
        double eps = sim_G_EPS(theSim);
        //double Sigma_nu = alpha*sim_GAMMALAW(theSim)*P*pow(r,1.5);
        //double Sigma_nu = alpha*sim_GAMMALAW(theSim)*P*(sqrt(M0)+sqrt(M1))/(sqrt(M0)*pow(dist_bh0,-1.5)+sqrt(M1)*pow(dist_bh1,-1.5));
        double Sigma_nu;
        if (sim_VISC_CONST(theSim)==1){
          Sigma_nu = rho*alpha;
        }else{
          Sigma_nu = alpha*P/sqrt(pow(dist_bh0*dist_bh0+eps*eps,-1.5)*M0+pow(dist_bh1*dist_bh1+eps*eps,-1.5)*M1);
        }
        c->cons[TAU] += dt*dV*Sigma_nu*rdrOm_a*(rdrOm_a + r*c->grad[UPP] + c->gradp[URR] );
        //c->cons[TAU] -= c->cons[TAU]*exp(-dist_bh0*dist_bh0/(.03*.03))/0.1*dt;
        //c->cons[TAU] -= c->cons[TAU]*exp(-dist_bh1*dist_bh1/(.03*.03))/0.1*dt;
      }
    }
  }
}

void cell_add_visc_src_old( struct Cell *** theCells ,struct Sim * theSim, double dt ){
  int i,j,k;
  for( k=0 ; k<sim_N(theSim,Z_DIR) ; ++k ){
    double zm = sim_FacePos(theSim,k-1,Z_DIR);
    double zp = sim_FacePos(theSim,k,Z_DIR);
    for( i=0 ; i<sim_N(theSim,R_DIR) ; ++i ){
      double rm = sim_FacePos(theSim,i-1,R_DIR);
      double rp = sim_FacePos(theSim,i,R_DIR);
      for( j=0 ; j<sim_N_p(theSim,i) ; ++j ){
        struct Cell * c = &(theCells[k][i][j]);
        double dphi = c->dphi;
        double rho = c->prim[RHO];
        double r   = .5*(rp+rm);
        double vr  = c->prim[URR];
        double dz = zp-zm;
        double dV = dphi*.5*(rp*rp-rm*rm)*dz;
        double nu = sim_EXPLICIT_VISCOSITY(theSim);
        c->cons[SRR] += -dt*dV*nu*rho*vr/r/r;
      }
    }
  }
}


