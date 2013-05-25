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

  double maxv = sqrt(cf2) + sqrt( vr*vr + vp*vp + vz*vz );

  return(maxv);

}

double cell_mindt( struct Cell *** theCells, struct Sim * theSim ){
  double dt_m = 1.e100;//HUGE_VAL;
  int i,j,k;
  for( k=sim_Nghost_min(theSim,Z_DIR) ; k<sim_N(theSim,Z_DIR)-sim_Nghost_max(theSim,Z_DIR) ; ++k ){
    double zm = sim_FacePos(theSim,k-1,Z_DIR);
    double zp = sim_FacePos(theSim,k,Z_DIR);
    double dz = zp-zm;
    for( i=sim_Nghost_min(theSim,R_DIR) ; i<sim_N(theSim,R_DIR)-sim_Nghost_max(theSim,R_DIR) ; ++i ){
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
        double Pp  = theCells[k][i][j].prim[PPP]; 


        double dt = sim_CFL(theSim)*dx/a;
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



