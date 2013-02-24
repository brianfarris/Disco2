#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdIO.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Grid.h"
#include "../Headers/Face.h"
#include "../Headers/GravMass.h"
#include "../Headers/header.h"


double fgrav( double M , double r , double eps ){
  double n = PHI_ORDER;
  return( M*pow(r,n-1.)/pow( pow(r,n) + pow(eps,n) ,1.+1./n) );
}

void gravMassaryForce( struct GravMass * theGravMasses , int p , double r , double phi , double * fr , double * fp ){

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

  double eps = G_EPS;

  double f1 = -fgrav( gravMass_M(theGravMasses,p) , script_r , eps );

  //double rH = theGravMasses[p].r*pow(theGravMasses[p].M/theGravMasses[0].M/3.,1./3.);
  //double pd = 0.8;
  //double fd = 1./(1.+exp(-( script_r/rH-pd)/(pd/10.)));

  //*fr = cosap*f1*fd;
  //*fp = sinap*f1*fd;
  *fr = cosap*f1;
  *fp = sinap*f1;
}

void source( double * prim , double * cons ,double divB, double * GradPsi, struct GravMass * theGravMasses , double rp , double rm , double phi , double dphi , double zp , double zm , double dt ){

  double rho = prim[RHO];
  double Pp  = prim[PPP];
  double r   = .5*(rp+rm);
  double vr  = prim[URR];
  double vz  = prim[UZZ];
  double vp  = prim[UPP]*r;
  double dz = zp-zm;
  double dV = dphi*.5*(rp*rp-rm*rm)*dz;

  double Br = prim[BRR];
  double Bp = prim[BPP];
  double Bz = prim[BZZ];

  double B2 = Br*Br+Bp*Bp+Bz*Bz;

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

  double vdotB = vr*Br+vp*Bp+vz*Bz;
  double BdotGradPsi = dV*(Br*GradPsi[0]+Bp*GradPsi[1]+Bz*GradPsi[2]);
  double vdotGradPsi = dV*(vr*GradPsi[0]+vp*GradPsi[1]+vz*GradPsi[2]);

  for( p=0 ; p<NP; ++p ){
    gravMassaryForce( theGravMasses , p , gravdist , phi , &fr , &fp );
    Fr += fr;
    Fp += fp;
  }
  cons[SRR] += dt*dV*( rho*vp*vp + Pp + .5*B2 - Bp*Bp )/r;
  cons[SRR] -= POWELL*dt*divB*Br;

  cons[SRR] += dt*dV*rho*Fr*sint;
  cons[SZZ] += dt*dV*rho*Fr*cost;
  cons[SZZ] -= POWELL*dt*divB*Bz;
  cons[LLL] += dt*dV*rho*Fp*r;
  cons[LLL] -= POWELL*dt*divB*Bp*r;
  cons[TAU] += dt*dV*rho*( Fr*(vr*sint+vz*cost) + Fp*vp );
  cons[TAU] -= POWELL*dt*(divB*vdotB+BdotGradPsi);

  double psi = prim[PSI];
  //cons[BRR] += dt*dV*( psi )/r;
  //cons[PSI] -= dt*dV*psi*DIVB_CH/DIVB_L;//DIVB_CH2/DIVB_CP2;
  cons[BRR] -= POWELL*dt*divB*vr/r;
  cons[BPP] -= POWELL*dt*divB*vp/r;
  cons[BZZ] -= POWELL*dt*divB*vz;
  cons[PSI] -= POWELL*dt*vdotGradPsi;

  if( INCLUDE_VISCOSITY ){
    double nu = EXPLICIT_VISCOSITY;
    cons[SRR] += -dt*dV*nu*rho*vr/r/r;
  }

}

void cell_add_src( struct Cell *** theCells ,struct Grid * theGrid, struct GravMass * theGravMasses , double dt ){
  int N_r_withghost = grid_N_r(theGrid)+grid_Nghost_rmin(theGrid)+grid_Nghost_rmax(theGrid);
  int N_z_withghost = grid_N_z(theGrid)+grid_Nghost_zmin(theGrid)+grid_Nghost_zmax(theGrid);

  int i,j,k;
  for( k=0 ; k<N_z_withghost ; ++k ){
    double zm = grid_z_faces(theGrid,k-1);
    double zp = grid_z_faces(theGrid,k);
    for( i=0 ; i<N_r_withghost ; ++i ){
      double rm = grid_r_faces(theGrid,i-1);
      double rp = grid_r_faces(theGrid,i);
      for( j=0 ; j<grid_N_p(theGrid,i) ; ++j ){
        struct Cell * c = &(theCells[k][i][j]);
        source( c->prim , c->cons , c->divB, c->GradPsi, theGravMasses , rp , rm , c->tiph-.5*c->dphi, c->dphi , zp , zm , dt );
      }
    }
  }
}


