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
  int i_m,j_m,k_m;
  double dt_m = 1.e100;//HUGE_VAL;
  double a_m,r_m;
  double mag_vel_m;
  int i,j,k;
  for( k=0 ; k<sim_N(theSim,Z_DIR) ; ++k ){
    double zm = sim_FacePos(theSim,k-1,Z_DIR);
    double zp = sim_FacePos(theSim,k,Z_DIR);
    double dz = zp-zm;
    for( i=0 ; i<sim_N(theSim,R_DIR) ; ++i ){
      double rm = sim_FacePos(theSim,i-1,R_DIR);
      double rp = sim_FacePos(theSim,i,R_DIR);
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
        double rho = theCells[k][i][j].prim[RHO]; 
        double Br  = theCells[k][i][j].prim[BRR];
        double Bp  = theCells[k][i][j].prim[BPP];
        double Bz  = theCells[k][i][j].prim[BZZ];
        double b2  = (Br*Br+Bp*Bp+Bz*Bz)/rho;

        double dt = sim_CFL(theSim)*dx/a;
        if( sim_EXPLICIT_VISCOSITY(theSim)>0.0 ){
          double dt_visc = .9*dx*dx/(sim_EXPLICIT_VISCOSITY(theSim));
          dt = dt/( 1. + dt/dt_visc );
        }
        if( dt_m > dt ) {
          mag_vel_m = sqrt(b2);
          a_m = a;
          dt_m = dt;
          i_m = i;
          j_m = j;
          k_m = k;
          r_m = r;
        }
      } 
    }
  }
  printf("r_m: %e, i_m: %d, j_m: %d, k_m: %d, a_m: %e mag_vel_m: %e\n",r_m,i_m,j_m,k_m,a_m,mag_vel_m);
  double dt2;
  MPI_Allreduce( &dt_m , &dt2 , 1 , MPI_DOUBLE , MPI_MIN , sim_comm );
  return( dt2 );

}



