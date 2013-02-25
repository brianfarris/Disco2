#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Grid.h"
#include "../Headers/Face.h"
#include "../Headers/GravMass.h"
#include "../Headers/header.h"



void cell_prim2cons(struct Cell ***theCells,struct Grid *theGrid) {
  int N_r_withghost = grid_N_r(theGrid)+grid_Nghost_rmin(theGrid)+grid_Nghost_rmax(theGrid);
  int N_z_withghost = grid_N_z(theGrid)+grid_Nghost_zmin(theGrid)+grid_Nghost_zmax(theGrid);

  int i,j,k;
  for(k = 0; k < N_z_withghost; k++){
    for(i = 0; i < N_r_withghost; i++){
      for(j = 0; j < grid_N_p(theGrid,i); j++){
        //printf("before r_plus: \n");
        double r_plus = grid_r_faces(theGrid,i);
        //printf("before r_minus:\n");
        double r_minus = grid_r_faces(theGrid,i-1);
        double r = 0.5*(r_minus+r_plus);
        //printf("before z_plus: \n");
        double z_plus = grid_z_faces(theGrid,k);
        //printf("before z_minus:\n");
        double z_minus = grid_z_faces(theGrid,k-1);
        double dz = z_plus-z_minus;
        double dphi = theCells[k][i][j].dphi;
        double dV = .5*(r_plus*r_plus - r_minus*r_minus)*dphi*dz;

        double rho = theCells[k][i][j].prim[RHO];
        double Pp  = theCells[k][i][j].prim[PPP];
        double vr  = theCells[k][i][j].prim[URR];
        double vp  = theCells[k][i][j].prim[UPP]*r;
        double vz  = theCells[k][i][j].prim[UZZ];

        double Br  = theCells[k][i][j].prim[BRR];
        double Bp  = theCells[k][i][j].prim[BPP];
        double Bz  = theCells[k][i][j].prim[BZZ];

        double v2  = vr*vr + vp*vp + vz*vz;
        double rhoe = Pp/(GAMMALAW - 1.);
        double B2 = Br*Br + Bp*Bp + Bz*Bz;

        theCells[k][i][j].cons[DDD] = rho*dV;
        theCells[k][i][j].cons[SRR] = rho*vr*dV;
        theCells[k][i][j].cons[LLL] = r*rho*vp*dV;
        theCells[k][i][j].cons[SZZ] = rho*vz*dV;
        theCells[k][i][j].cons[TAU] = ( .5*rho*v2 + rhoe + .5*B2 )*dV;

        theCells[k][i][j].cons[BRR] = Br*dV/r;
        theCells[k][i][j].cons[BPP] = Bp*dV/r;
        theCells[k][i][j].cons[BZZ] = Bz*dV;

        theCells[k][i][j].cons[PSI] = theCells[k][i][j].prim[PSI]*dV;
      }
    }
  }
}


void cons2prim( double * cons , double * prim , double r , double dV ){

  double rho = cons[DDD]/dV;
  if( rho < RHO_FLOOR ) rho = RHO_FLOOR;
  double Sr  = cons[SRR]/dV;
  double Sp  = cons[LLL]/dV/r;
  double Sz  = cons[SZZ]/dV;
  double E   = cons[TAU]/dV;

  double Br  = cons[BRR]/dV*r;
  double Bp  = cons[BPP]/dV*r;
  double Bz  = cons[BZZ]/dV;

  double B2 = Br*Br+Bp*Bp+Bz*Bz;
  double vr = Sr/rho;
  double vp = Sp/rho;
  double vz = Sz/rho;
  double KE = .5*( Sr*vr + Sp*vp + Sz*vz );
  double rhoe = E - KE - .5*B2;
  double Pp = (GAMMALAW-1.)*rhoe;
  double v_magnitude = sqrt(2.*KE/rho);

  if( Pp < CS_FLOOR*CS_FLOOR*rho/GAMMALAW ) Pp = CS_FLOOR*CS_FLOOR*rho/GAMMALAW;
  if ( Pp > CS_CAP*CS_CAP*rho/GAMMALAW) Pp = CS_CAP*CS_CAP*rho/GAMMALAW;
  if ( v_magnitude>VEL_CAP ){
    vr = vr*VEL_CAP/v_magnitude;
    vp = vp*VEL_CAP/v_magnitude;
    vz = vz*VEL_CAP/v_magnitude;
  }
  prim[RHO] = rho;
  prim[PPP] = Pp;
  prim[URR] = vr;
  prim[UPP] = vp/r;
  prim[UZZ] = vz;

  prim[BRR] = Br;
  prim[BPP] = Bp;
  prim[BZZ] = Bz;

  prim[PSI] = cons[PSI]/dV;

  /*
     int q;
     for( q=NUM_C ; q<NUM_Q ; ++q ){
     prim[q] = cons[q]/cons[DDD];
     }
     */
}


void cell_calc_prim( struct Cell *** theCells ,struct Grid * theGrid){
  int N_r_withghost = grid_N_r(theGrid)+grid_Nghost_rmin(theGrid)+grid_Nghost_rmax(theGrid);
  int N_z_withghost = grid_N_z(theGrid)+grid_Nghost_zmin(theGrid)+grid_Nghost_zmax(theGrid);

  int i,j,k;
  for( k=0 ; k<N_z_withghost ; ++k ){
    double zm = grid_z_faces(theGrid,k-1);
    double zp = grid_z_faces(theGrid,k);
    double dz = zp-zm;
    for( i=0 ; i<N_r_withghost ; ++i ){
      double rm = grid_r_faces(theGrid,i-1);
      double rp = grid_r_faces(theGrid,i);
      double r = .5*(rp+rm);
      for( j=0 ; j<grid_N_p(theGrid,i) ; ++j ){
        struct Cell * c = &(theCells[k][i][j]);
        double dV = .5*(rp*rp-rm*rm)*c->dphi*dz;
        cons2prim( c->cons , c->prim , r , dV );
      }
    }
  }
}



