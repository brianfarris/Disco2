#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Grid.h"
#include "../Headers/header.h"

double maxvel(double * prim , double w , double r ,struct Grid * theGrid){

  double Pp  = prim[PPP];
  double rho = prim[RHO];
  double vp  = prim[UPP]*r-w;
  double vr  = prim[URR];
  double vz  = prim[UZZ];
  double cs  = sqrt(grid_GAMMALAW(theGrid)*Pp/rho);
  double cf2 = cs*cs;

  if (grid_runtype(theGrid)==MHD){
    double Br  = prim[BRR];
    double Bp  = prim[BPP];
    double Bz  = prim[BZZ];
    double b2  = (Br*Br+Bp*Bp+Bz*Bz)/rho;
    cf2 += b2;
  }
  double maxv = sqrt(cf2) + sqrt( vr*vr + vp*vp + vz*vz );

  if (grid_runtype(theGrid)==MHD){
    double ch = grid_DIVB_CH(theGrid); 
    if( maxv < ch ) maxv = ch;
  }

  return(maxv);

}

double cell_mindt( struct Cell *** theCells, struct Grid * theGrid ){
  int N_r_withghost = grid_N_r(theGrid)+grid_Nghost_rmin(theGrid)+grid_Nghost_rmax(theGrid);
  int N_z_withghost = grid_N_z(theGrid)+grid_Nghost_zmin(theGrid)+grid_Nghost_zmax(theGrid);

  double dt_m = HUGE_VAL;
  int i,j,k;
  for( k=0 ; k<N_z_withghost ; ++k ){
    double zm = grid_z_faces(theGrid,k-1);
    double zp = grid_z_faces(theGrid,k);
    double dz = zp-zm;
    for( i=0 ; i<N_r_withghost ; ++i ){
      double rm = grid_r_faces(theGrid,i-1);
      double rp = grid_r_faces(theGrid,i);
      double dr = rp-rm;
      double r = .5*(rp+rm);
      for( j=0 ; j<grid_N_p(theGrid,i) ; ++j ){
        int jm = j-1;
        if( j==0 ) jm = grid_N_p(theGrid,i)-1;
        double w = .5*(theCells[k][i][j].wiph+theCells[k][i][jm].wiph);
        double dx = dr;
        double rdphi = .5*(rp+rm)*theCells[k][i][j].dphi;
        if( rdphi<dr ) dx = rdphi;
        if( dx>dz ) dx = dz;
        double a  = maxvel( theCells[k][i][j].prim , w , r ,theGrid);
        double dt = grid_CFL(theGrid)*dx/a;
        if( grid_INCLUDE_VISCOSITY(theGrid)==1 ){
          double dt_visc = .9*dx*dx/(grid_EXPLICIT_VISCOSITY(theGrid));
          dt = dt/( 1. + dt/dt_visc );
        }
        if( dt_m > dt ) dt_m = dt;
      } 
    }
  }
  double dt2;
  MPI_Allreduce( &dt_m , &dt2 , 1 , MPI_DOUBLE , MPI_MIN , grid_comm );
  return( dt2 );

}



