#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Sim.h"
#include "../Headers/header.h"

double maxvel(double * prim , double w , double r ,struct Sim * theSim){

  double Pp  = prim[PPP];
  double rho = prim[RHO];
  double vp  = prim[UPP]*r-w;
  double vr  = prim[URR];
  double vz  = prim[UZZ];
  double cs  = sqrt(sim_GAMMALAW(theSim)*Pp/rho);
  double cf2 = cs*cs;

  if (sim_runtype(theSim)==MHD){
    double Br  = prim[BRR];
    double Bp  = prim[BPP];
    double Bz  = prim[BZZ];
    double b2  = (Br*Br+Bp*Bp+Bz*Bz)/rho;
    cf2 += b2;
  }
  double maxv = sqrt(cf2) + sqrt( vr*vr + vp*vp + vz*vz );

  if (sim_runtype(theSim)==MHD){
    double ch = sim_DIVB_CH(theSim); 
    if( maxv < ch ) maxv = ch;
  }

  return(maxv);

}

double cell_mindt( struct Cell *** theCells, struct Sim * theSim ){

  double dt_m = 1.e100;//HUGE_VAL;
  int i,j,k;
  for( k=0 ; k<sim_N_z(theSim) ; ++k ){
    double zm = sim_z_faces(theSim,k-1);
    double zp = sim_z_faces(theSim,k);
    double dz = zp-zm;
    for( i=0 ; i<sim_N_r(theSim) ; ++i ){
      double rm = sim_r_faces(theSim,i-1);
      double rp = sim_r_faces(theSim,i);
      double dr = rp-rm;
      double r = .5*(rp+rm);
      for( j=0 ; j<sim_N_p(theSim,i) ; ++j ){
        int jm = j-1;
        if( j==0 ) jm = sim_N_p(theSim,i)-1;
        double w = .5*(theCells[k][i][j].wiph+theCells[k][i][jm].wiph);
        double dx = dr;
        double rdphi = .5*(rp+rm)*theCells[k][i][j].dphi;
        if( rdphi<dr ) {
          dx = rdphi;
        }
        if( dx>dz ) {
          dx = dz;
        }
        double a  = maxvel( theCells[k][i][j].prim , w , r ,theSim);
        double dt = sim_CFL(theSim)*dx/a;
        if( sim_INCLUDE_VISCOSITY(theSim)==1 ){
          double dt_visc = .9*dx*dx/(sim_EXPLICIT_VISCOSITY(theSim));
          dt = dt/( 1. + dt/dt_visc );
        }
        if( dt_m > dt ) {
          dt_m = dt;
        }
      } 
    }
  }
  double dt2;
  MPI_Allreduce( &dt_m , &dt2 , 1 , MPI_DOUBLE , MPI_MIN , sim_comm );
  return( dt2 );

}



