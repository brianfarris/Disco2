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
 
  double one_o_nu = 1.0/(alpha*P/rho)*sqrt(pow(r0,-3)*M0/Mtotal + pow(r1,-3)*M1/Mtotal);
  double t_visc0 = TVISC_FAC * 2./3. * r0*r0*one_o_nu;
  double t_visc1 = TVISC_FAC * 2./3. * r1*r1*one_o_nu;    
 
  if (t_visc0 < 10.* dt){
    t_visc0 = 10.*dt;
  }
  if (t_visc1 < 10.* dt){
    t_visc1 = 10.*dt;
  }

  *drho_dt_sink = 0.0;
  double sink_size = 0.5;

  if (p==0){
    if (r0<sink_size){
      *drho_dt_sink = rho / t_visc0;
    }
  } else if (p==1){
    if (r1<sink_size){
      *drho_dt_sink = rho / t_visc1;
    }
  } else{
    printf("bad value for p\n");
    exit(1);
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

  double total_torque_temp = 0.0;
  //  double total_Fr_temp = 0.0;


  double Mtotal = 1.0;
  double sep = 1.0;
  double M0 = gravMass_M(theGravMasses,0);
  double M1 = gravMass_M(theGravMasses,1);

  double rbh0 = gravMass_r(theGravMasses,0);
  double rbh1 = gravMass_r(theGravMasses,1);

  double phi_bh0 = gravMass_phi(theGravMasses,0);
  double phi_bh1 = gravMass_phi(theGravMasses,1);

  double xbh0 = rbh0*cos(phi_bh0);
  double ybh0 = rbh0*sin(phi_bh0);

  double xbh1 = rbh1*cos(phi_bh1);
  double ybh1 = rbh1*sin(phi_bh1);



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
        double dV = dphi*.5*(rp*rp-rm*rm)*dz;

        double xpos = r*cos(phi);
        double ypos = r*sin(phi);

        double r0 = sqrt((xpos-xbh0)*(xpos-xbh0)+(ypos-ybh0)*(ypos-ybh0));
        double r1 = sqrt((xpos-xbh1)*(xpos-xbh1)+(ypos-ybh1)*(ypos-ybh1));

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
        int p;
        for( p=0 ; p<sim_NumGravMass(theSim); ++p ){
          gravMassForce( theGravMasses , theSim , p , gravdist , phi , &fr , &fp );
          Fr += fr;
          Fp += fp;
        }
        double w_a = sim_W_A(theSim,r);
        double drOm_a = sim_OM_A_DERIV(theSim,r);
        double F_centrifugal_r = w_a*w_a/r;
        double F_coriolis_r =  2.*w_a*vp/r;
        double F_coriolis_phi = -2.*w_a*vr/r;
        double F_euler_phi =  -vr*r*drOm_a;

        c->cons[SRR] += dt*dV*( rho*vp*vp + Pp )/r;
        c->cons[SRR] += dt*dV*rho*(Fr*sint+F_centrifugal_r/*+F_coriolis_r*/);
        c->cons[LLL] += dt*dV*rho*(Fp/*+F_coriolis_phi+F_euler_phi*/)*r;
        c->cons[SZZ] += dt*dV*rho*Fr*cost;
        c->cons[TAU] += dt*dV*rho*( (Fr*sint+F_centrifugal_r)*vr+Fr*vz*cost + (Fp+F_euler_phi)*vp);
        if ((i>=imin_noghost) && (i<imax_noghost) && (k>=kmin_noghost) && (k<kmax_noghost)){
          total_torque_temp += Fp*r;
          //total_Fr_temp += Fr;
        }

        if (sim_RhoSinkOn(theSim)==1){
          double drho_dt_sink;
          for( p=0 ; p<sim_NumGravMass(theSim); ++p ){
            get_rho_sink( theGravMasses, theSim,p,dt,r0,r1,M0,M1,rho,Pp, &drho_dt_sink);
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
          double vdotGradPsi = /*dV**/(vr*GradPsi[0]+(vp-sim_W_A(theSim,r))*GradPsi[1]+vz*GradPsi[2]);

          c->cons[SRR] += dt*dV*( .5*B2 - Bp*Bp )/r;
          c->cons[TAU] += dt*dV*r*Br*Bp*drOm_a;
          double psi = c->prim[PSI];
          c->cons[BRR] += dt*dV*( psi )/r;
          //cons[PSI] -= dt*dV*psi*DIVB_CH/DIVB_L;//DIVB_CH2/DIVB_CP2;
          //c->cons[BPP] += dt*dV*Br*sim_OM_A_DERIV(theSim,r);
          c->cons[BPP] += dt*dV*Br*drOm_a;

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
          double cooling_cap = c->cons[TAU]*1000000000000000.0;
          //double cooling_term = numfac * pow(Pp/rho,0.5)*pow(r,-1.5)*pow(a_o_M,0.5)*dt*dV;
          double cooling_term = 9./4.*0.1*5./3.*1.e6/rho*pow(Pp/rho,4.)*dt*dV;
          if (cooling_term>cooling_cap) cooling_term = cooling_cap;
          c->cons[TAU] -= cooling_term;
        }
      }
    }
  }

  double Mdot_reduce[2];
  MPI_Allreduce( Mdot_temp,Mdot_reduce , 2, MPI_DOUBLE, MPI_SUM, sim_comm);
  double total_torque_reduce;
  MPI_Allreduce( &total_torque_temp,&total_torque_reduce, 1, MPI_DOUBLE, MPI_SUM, sim_comm);
  int p;
  for( p=0 ; p<sim_NumGravMass(theSim); ++p ){
    gravMass_set_Mdot(theGravMasses,Mdot_reduce[p],p);
    gravMass_set_total_torque(theGravMasses,total_torque_reduce,p);
  }
}

void cell_add_visc_src( struct Cell *** theCells ,struct Sim * theSim, struct GravMass * theGravMasses, double dt ){
  double Mtotal = 1.0;
  double sep = 1.0;
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
        double drOm_a = sim_OM_A_DERIV(theSim,r);
        //double Sigma_nu = alpha*sim_GAMMALAW(theSim)*P*pow(r,1.5);
        //double Sigma_nu = alpha*sim_GAMMALAW(theSim)*P*(sqrt(M0)+sqrt(M1))/(sqrt(M0)*pow(dist_bh0,-1.5)+sqrt(M1)*pow(dist_bh1,-1.5));
        double Sigma_nu = alpha*P/sqrt(pow(dist_bh0,-3)*M0+pow(dist_bh1,-3.)*M1);

        c->cons[TAU] += dt*dV*Sigma_nu*r*drOm_a*(r*drOm_a + r*c->grad[UPP] + c->gradp[URR]);
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


