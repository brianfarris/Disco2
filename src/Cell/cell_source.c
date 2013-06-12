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
  double Om = 20.0;
  return( M*r*Om*Om );
}

void get_rho_sink( struct GravMass * theGravMasses, double RhoSinkTimescale, int p, double r, double phi, double rho, double * rho_sink){
  double rp = gravMass_r(theGravMasses,p);
  double pp = gravMass_phi(theGravMasses,p);
  double cosp = cos(phi);
  double sinp = sin(phi);
  double dx = r*cosp-rp*cos(pp);
  double dy = r*sinp-rp*sin(pp);
  double script_r = sqrt(dx*dx+dy*dy);

  *rho_sink = rho / (RhoSinkTimescale*2.0*M_PI) * exp(-script_r*script_r/(0.25*0.25));

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
  //FILE * gradpsifile= fopen("gradpsi.dat","w");
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
        c->cons[SRR] += dt*dV*( rho*vp*vp + Pp )/r;
        c->cons[SRR] += dt*dV*rho*Fr*sint;
        c->cons[LLL] += dt*dV*rho*Fp*r;
        c->cons[SZZ] += dt*dV*rho*Fr*cost;
        c->cons[TAU] += dt*dV*rho*( Fr*(vr*sint+vz*cost) + Fp*vp );

        if (sim_RhoSinkTimescale(theSim)>0.0){
          double rho_sink;
          for( p=0 ; p<sim_NumGravMass(theSim); ++p ){
            get_rho_sink(theGravMasses,sim_RhoSinkTimescale(theSim),p,gravdist,phi,rho, &rho_sink);
            c->cons[RHO] -= rho_sink * dt * dV;
          }
        }

        if(sim_runtype(theSim)==MHD){
          double Br = c->prim[BRR];
          double Bp = c->prim[BPP];
          double Bz = c->prim[BZZ];
          double B2 = Br*Br+Bp*Bp+Bz*Bz;
          double divB = c->divB;
          double GradPsi[3];
          //GradPsi[0] = 0.0;//r*c->GradPsi[0];
          //GradPsi[1] = 0.0;//r*c->GradPsi[1];
          //GradPsi[2] = 0.0;//r*c->GradPsi[2];
          GradPsi[0] = r*c->GradPsi[0];
          GradPsi[1] = r*c->GradPsi[1];
          GradPsi[2] = r*c->GradPsi[2];

          double vdotB = vr*Br+vp*Bp+vz*Bz;
          double BdotGradPsi = /*dV**/(Br*GradPsi[0]+Bp*GradPsi[1]+Bz*GradPsi[2]);
          double vdotGradPsi = /*dV**/(vr*GradPsi[0]+(vp-sim_W_A(theSim,r))*GradPsi[1]+vz*GradPsi[2]);

          /*
             if (fabs(z)<0.00001){
             fprintf(gradpsifile,"%e %e %e %e %e\n",r,phi,dV,GradPsi[0]/dV, c->prim[PSI]);
             }
             */
          c->cons[SRR] += dt*dV*( .5*B2 - Bp*Bp )/r;
          c->cons[SRR] -= POWELL*dt*divB*Br;
          c->cons[SZZ] -= POWELL*dt*divB*Bz;
          c->cons[LLL] -= POWELL*dt*divB*Bp*r;
          c->cons[TAU] -= POWELL*dt*(divB*vdotB+BdotGradPsi);

          //double psi = c->prim[PSI];
          //cons[BRR] += dt*dV*( psi )/r;
          //cons[PSI] -= dt*dV*psi*DIVB_CH/DIVB_L;//DIVB_CH2/DIVB_CP2;
          c->cons[BRR] -= POWELL*dt*divB*vr/r;
          c->cons[BPP] -= POWELL*dt*divB*(vp-sim_W_A(theSim,r))/r;
          c->cons[BZZ] -= POWELL*dt*divB*vz;
          c->cons[PSI] -= POWELL*dt*vdotGradPsi;

          c->cons[BPP] += dt*dV*Br*sim_OM_A_DERIV(theSim,r);

        }
      }
    }
  }
  //fclose(gradpsifile);

}

void cell_add_visc_src( struct Cell *** theCells ,struct Sim * theSim, double dt ){
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
        // c->cons[SRR] += -dt*dV*nu*rho*vr/r/r;
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


