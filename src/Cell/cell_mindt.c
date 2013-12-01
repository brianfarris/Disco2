#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Metric.h"
#include "../Headers/Sim.h"
#include "../Headers/header.h"

double maxvel(double * prim , double w , double r ,struct Sim * theSim){

    double maxv;

    double Pp  = prim[PPP];
    double rho = prim[RHO];
    double vp  = prim[UPP]*r-w;
    double vr  = prim[URR];
    double vz  = prim[UZZ];
    double cf2 = sim_GAMMALAW(theSim)*Pp/rho;

    maxv = sqrt(cf2) + sqrt( vr*vr + vp*vp + vz*vz );

    return maxv;
}

double cell_maxvel_gr(double *prim, int dir, double w, double r, struct Sim *theSim)
{
    double maxv;
    double rho, Pp, v[3], GAMMALAW, rhoh, vn;
    double a, bn, b[3], u02;
    double vp, vm, sig, dv, gam, cs2;
    int i;
    struct Metric *g;

    rho = prim[RHO];
    Pp  = prim[PPP];
    v[0] = prim[URR];
    v[1] = prim[UPP];
    v[2] = prim[UZZ];
    vn = v[dir];

    GAMMALAW = sim_GAMMALAW(theSim);
    rhoh = rho + GAMMALAW * Pp / (GAMMALAW-1.0);
    
    g = metric_create(time_global, r, 0, 0);
    a = metric_lapse(g);
    for(i=0; i<3; i++)
        b[i] = metric_shift_u(g,i);
    bn = b[dir];
    gam = metric_gamma_uu(g,dir,dir);
    u02 = 1.0/(-metric_g_dd(g,0,0)-2.0*metric_dot3_u(g,b,v)-metric_square3_u(g,v));

    cs2 = GAMMALAW*Pp/rhoh;
    sig = cs2 / ((1-cs2)*u02*a*a);
    dv = sqrt(sig*(1.0+sig)*a*a*gam - sig*(vn+bn)*(vn+bn));

    if(dir == PDIRECTION)
    {
        vm = fabs((vn - sig*bn - dv) / (1.0+sig) - w);
        vp = fabs((vn - sig*bn + dv) / (1.0+sig) - w);
    }
    else
    {
        vm = fabs((vn - sig*bn - dv) / (1.0+sig));
        vp = fabs((vn - sig*bn + dv) / (1.0+sig));
    }

    if(vm < vp)
        maxv = vp;
    else
        maxv = vm;

    metric_destroy(g);

    return maxv;
}

double cell_mindt_newt( struct Cell *** theCells, struct Sim * theSim ){
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

double cell_mindt_gr(struct Cell ***theCells, struct Sim *theSim)
{
    double dt_m = 1.e100;//HUGE_VAL;
    int i,j,k;
    for(k = sim_Nghost_min(theSim,Z_DIR); k < sim_N(theSim,Z_DIR)-sim_Nghost_max(theSim,Z_DIR); ++k)
    {
        double zm = sim_FacePos(theSim,k-1,Z_DIR);
        double zp = sim_FacePos(theSim,k,Z_DIR);
        double dz = zp-zm;
        
        for(i = sim_Nghost_min(theSim,R_DIR); i < sim_N(theSim,R_DIR)-sim_Nghost_max(theSim,R_DIR); ++i)
        {
            double rm = sim_FacePos(theSim,i-1,R_DIR);
            double rp = sim_FacePos(theSim,i,R_DIR);
            double dr = rp-rm;
            double r = .5*(rp+rm);
            
            for(j = 0; j < sim_N_p(theSim,i); ++j)
            {
                int jm = j-1;
                if(j == 0) 
                    jm = sim_N_p(theSim,i)-1;
                double w = .5*(theCells[k][i][j].wiph+theCells[k][i][jm].wiph); 
                double dphi = theCells[k][i][j].dphi;

                double ar, ap, az;
                ar = cell_maxvel_gr(theCells[k][i][j].prim, 0, w, r, theSim);
                ap = cell_maxvel_gr(theCells[k][i][j].prim, 1, w, r, theSim);
                az = cell_maxvel_gr(theCells[k][i][j].prim, 2, w, r, theSim);
                
                double dt = dr/ar;
                if(dt > dphi/ap)
                    dt = dphi/ap;
                if(dt > dz/az)
                    dt = dz/az;
                dt *= sim_CFL(theSim);

                if(dt_m > dt)
                    dt_m = dt;
            }     
        }
    }
    double dt2;
    MPI_Allreduce( &dt_m , &dt2 , 1 , MPI_DOUBLE , MPI_MIN , sim_comm );

    double dt_newt = cell_mindt_newt(theCells, theSim);
    //printf("\nNewtonian dt: %lg, GR dt: %lg\n", dt_newt, dt2);

    return dt2;
}

