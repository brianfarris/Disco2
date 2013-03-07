#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Grid.h"
#include "../Headers/Face.h"
#include "../Headers/GravMass.h"
#include "../Headers/header.h"


double fgrav( double M , double r , double eps, double n ){
  return( M*pow(r,n-1.)/pow( pow(r,n) + pow(eps,n) ,1.+1./n) );
}

void gravMassForce( struct GravMass * theGravMasses ,struct Grid * theGrid, int p , double r , double phi , double * fr , double * fp ){
  double G_EPS=grid_G_EPS(theGrid);
  double PHI_ORDER=grid_PHI_ORDER(theGrid);

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

  double f1 = -fgrav( gravMass_M(theGravMasses,p) , script_r , G_EPS, PHI_ORDER );

  //double rH = theGravMasses[p].r*pow(theGravMasses[p].M/theGravMasses[0].M/3.,1./3.);
  //double pd = 0.8;
  //double fd = 1./(1.+exp(-( script_r/rH-pd)/(pd/10.)));

  //*fr = cosap*f1*fd;
  //*fp = sinap*f1*fd;
  *fr = cosap*f1;
  *fp = sinap*f1;
}

void cell_add_src( struct Cell *** theCells ,struct Grid * theGrid, struct GravMass * theGravMasses , double dt ){
  int GRAV2D=grid_GRAV2D(theGrid);
  int POWELL=grid_POWELL(theGrid);
  int i,j,k;
  for( k=0 ; k<grid_N_z(theGrid) ; ++k ){
    double zm = grid_z_faces(theGrid,k-1);
    double zp = grid_z_faces(theGrid,k);
    for( i=0 ; i<grid_N_r(theGrid) ; ++i ){
      double rm = grid_r_faces(theGrid,i-1);
      double rp = grid_r_faces(theGrid,i);
      for( j=0 ; j<grid_N_p(theGrid,i) ; ++j ){
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
        for( p=0 ; p<grid_NumGravMass(theGrid); ++p ){
          gravMassForce( theGravMasses , theGrid , p , gravdist , phi , &fr , &fp );
          Fr += fr;
          Fp += fp;
        }
        c->cons[SRR] += dt*dV*( rho*vp*vp + Pp )/r;
        c->cons[SRR] += dt*dV*rho*Fr*sint;
        c->cons[LLL] += dt*dV*rho*Fp*r;
        c->cons[SZZ] += dt*dV*rho*Fr*cost;
        c->cons[TAU] += dt*dV*rho*( Fr*(vr*sint+vz*cost) + Fp*vp );

        if(grid_runtype(theGrid)==MHD){
          double Br = c->prim[BRR];
          double Bp = c->prim[BPP];
          double Bz = c->prim[BZZ];
          double B2 = Br*Br+Bp*Bp+Bz*Bz;
          double divB = c->divB;
          double *GradPsi = c->GradPsi;
          double vdotB = vr*Br+vp*Bp+vz*Bz;
          double BdotGradPsi = dV*(Br*GradPsi[0]+Bp*GradPsi[1]+Bz*GradPsi[2]);
          double vdotGradPsi = dV*(vr*GradPsi[0]+vp*GradPsi[1]+vz*GradPsi[2]);

          c->cons[SRR] += dt*dV*( .5*B2 - Bp*Bp )/r;
          c->cons[SRR] -= POWELL*dt*divB*Br;
          c->cons[SZZ] -= POWELL*dt*divB*Bz;
          c->cons[LLL] -= POWELL*dt*divB*Bp*r;
          c->cons[TAU] -= POWELL*dt*(divB*vdotB+BdotGradPsi);

          double psi = c->prim[PSI];
          //cons[BRR] += dt*dV*( psi )/r;
          //cons[PSI] -= dt*dV*psi*DIVB_CH/DIVB_L;//DIVB_CH2/DIVB_CP2;
          c->cons[BRR] -= POWELL*dt*divB*vr/r;
          c->cons[BPP] -= POWELL*dt*divB*vp/r;
          c->cons[BZZ] -= POWELL*dt*divB*vz;
          c->cons[PSI] -= POWELL*dt*vdotGradPsi;
        }
      }
    }
  }
}

void cell_add_visc_src( struct Cell *** theCells ,struct Grid * theGrid, double dt ){
  int INCLUDE_VISCOSITY=grid_INCLUDE_VISCOSITY(theGrid);
  double EXPLICIT_VISCOSITY=grid_EXPLICIT_VISCOSITY(theGrid);
  int i,j,k;
  for( k=0 ; k<grid_N_z(theGrid) ; ++k ){
    double zm = grid_z_faces(theGrid,k-1);
    double zp = grid_z_faces(theGrid,k);
    for( i=0 ; i<grid_N_r(theGrid) ; ++i ){
      double rm = grid_r_faces(theGrid,i-1);
      double rp = grid_r_faces(theGrid,i);
      for( j=0 ; j<grid_N_p(theGrid,i) ; ++j ){
        struct Cell * c = &(theCells[k][i][j]);
        double dphi = c->dphi;
        double rho = c->prim[RHO];
        double r   = .5*(rp+rm);
        double vr  = c->prim[URR];
        double dz = zp-zm;
        double dV = dphi*.5*(rp*rp-rm*rm)*dz;
         if (grid_INCLUDE_VISCOSITY(theGrid)){
          double nu = grid_EXPLICIT_VISCOSITY(theGrid);
          c->cons[SRR] += -dt*dV*nu*rho*vr/r/r;
        }
      }
    }
  }
}


