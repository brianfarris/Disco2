#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/EOS.h"
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

double cell_maxvel_gr(double *prim, int dir, double w, double r, double phi, double z, struct Sim *theSim)
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
    
    g = metric_create(time_global, r, phi, z, theSim);
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

//TODO: WRITE THIS.
double cell_maxvel_grdisc(double *prim, int dir, double w, double r, double phi, double z, struct Sim *theSim)
{
    return 0.0;
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
    double M = sim_GravM(theSim);
    int i,j,k;
    for(k = sim_Nghost_min(theSim,Z_DIR); k < sim_N(theSim,Z_DIR)-sim_Nghost_max(theSim,Z_DIR); ++k)
    {
        double zm = sim_FacePos(theSim,k-1,Z_DIR);
        double zp = sim_FacePos(theSim,k,Z_DIR);
        double dz = zp-zm;
        double z = 0.5*(zm+zp);
        
        for(i = sim_Nghost_min(theSim,R_DIR); i < sim_N(theSim,R_DIR)-sim_Nghost_max(theSim,R_DIR); ++i)
        {
            double rm = sim_FacePos(theSim,i-1,R_DIR);
            double rp = sim_FacePos(theSim,i,R_DIR);
            double dr = rp-rm;
            double r = .5*(rp+rm);
            
            //TODO: cleaner way of doing this?
            if(r < 2*M)
                continue;

            for(j = 0; j < sim_N_p(theSim,i); ++j)
            {
                int jm = j-1;
                if(j == 0) 
                    jm = sim_N_p(theSim,i)-1;
                double w = .5*(theCells[k][i][j].wiph+theCells[k][i][jm].wiph); 
                double dphi = theCells[k][i][j].dphi;
                double phi = theCells[k][i][j].tiph - 0.5*dphi;

                double ar, ap, az;
                ar = cell_maxvel_gr(theCells[k][i][j].prim, 0, w, r, phi, z, theSim);
                ap = cell_maxvel_gr(theCells[k][i][j].prim, 1, w, r, phi, z, theSim);
                az = cell_maxvel_gr(theCells[k][i][j].prim, 2, w, r, phi, z, theSim);
                
                double dx = dr;
                double dt = dr/ar;
                if(dt > dphi/ap)
                {
                    dt = dphi/ap;
                    dx = r*dphi;
                }
                if(dt > dz/az)
                {
                    dt = dz/az;
                    dx = dz;
                }
                dt *= sim_CFL(theSim);

                if(sim_Background(theSim) == GRVISC1)
                {
                    //TODO: Incorporate proper alpha viscosity
                    double rho = theCells[k][i][j].prim[RHO];
                    double Pp = theCells[k][i][j].prim[PPP];
                    double GAM = sim_GAMMALAW(theSim);
                    double rhoh = rho + GAM/(GAM-1)*Pp;
                    double nu;
                    double alpha = sim_AlphaVisc(theSim);
                    double height = 2.0 * sqrt(r*r*r*(1-3.0*M/r)/M) * sqrt(GAM*Pp/rhoh);
                    if(alpha < 0.0)
                        nu = -alpha;
                    else
                    {
                        nu = alpha * sqrt(GAM*Pp/rhoh) * height;
                    }
                    double dtv = 0.25*dx*dx/nu;
                    
                    //Luminosity Time
                    double maxdisp = 0.1; //Should be a fraction like 0.1
                    double dtl;

                    if(sim_CoolingType(theSim) != COOL_NONE)
                    {
                        double dV = 0.5*(rp+rm)*dr*dphi*dz;
                        double rho = theCells[k][i][j].prim[RHO];
                        double Pp = theCells[k][i][j].prim[PPP];
                        double U[4];
                        double cool;

                        cool = eos_cool(theCells[k][i][j].prim, height, theSim);
                    /*
                        if(sim_CoolingType(theSim) == COOL_ISOTHERM)
                        {
                            cool = (Pp/rho - sim_CoolPar1(theSim)) / sim_CoolPar2(theSim);
                        }
                        else if(sim_CoolingType(theSim) == COOL_BB_ES)
                        {
                            double temp = Pp/rho;
                            cool = sim_CoolPar1(theSim) * pow(temp, 4) / rho;
                        }
                    */
                        struct Metric *g;
                        int mu;
                        g = metric_create(time_global, r, phi, z, theSim);
                        double al, b[3], sqrtg, v[3], u[4], u_d[4];
                        al = metric_lapse(g);
                        for(mu=0; mu<3; mu++)
                            b[mu] = metric_shift_u(g,mu);
                        sqrtg = metric_sqrtgamma(g)/r;
                        for(mu = 0; mu < 4; mu++)
                            U[mu] = metric_frame_U_u(g, mu, theSim);
                        v[0] = theCells[k][i][j].prim[URR];
                        v[1] = theCells[k][i][j].prim[UPP];
                        v[2] = theCells[k][i][j].prim[UZZ];
                        u[0] = 1.0/sqrt(-metric_g_dd(g,0,0) - 2*metric_dot3_u(g,b,v) - metric_square3_u(g,v));
                        u[1] = u[0]*v[0]; u[2] = u[0]*v[1]; u[3] = u[0]*v[2];
                        for(mu=0; mu<4; mu++)
                            u_d[mu] = 0.0;
                        for(mu=0; mu<4; mu++)
                        {
                            u_d[0] += metric_g_dd(g,0,mu)*u[mu];
                            u_d[1] += metric_g_dd(g,1,mu)*u[mu];
                            u_d[2] += metric_g_dd(g,2,mu)*u[mu];
                            u_d[3] += metric_g_dd(g,3,mu)*u[mu];
                        }
                        metric_destroy(g);

                        double dtl_min;
                        dtl = maxdisp * fabs(theCells[k][i][j].cons[TAU] / (dV*sqrtg*al* cool * (u_d[0]*U[0]+u_d[1]*U[1]+u_d[2]*U[2]+u_d[3]*U[3])));
                        dtl = fabs(dtl);
                        if (dtl < 0)
                            printf("WHAT.");
                        dtl_min = dtl;
                        dtl = maxdisp * fabs(theCells[k][i][j].cons[SRR] / (dV*al*sqrtg*u_d[1]*cool));
                        if(dtl < dtl_min)
                            dtl_min = dtl;
                        dtl = maxdisp * fabs(theCells[k][i][j].cons[LLL] / (dV*al*sqrtg*u_d[2]*cool));
                        if(dtl < dtl_min)
                            dtl_min = dtl;
                        if(sim_N(theSim, Z_DIR) > 1)
                        {
                            dtl = maxdisp * fabs(theCells[k][i][j].cons[SZZ] / (dV*al*sqrtg*u_d[3]*cool));
                            if(dtl < dtl_min)
                                dtl_min = dtl;
                        }
                        dtl = dtl_min;

                        dt = dt/(1.0 + dt/dtv + dt/dtl);
                    }
                    else
                    {
                        dt = dt/(1.0 + dt/dtv);
                    }

                }

                if(dt_m > dt)
                    dt_m = dt;
            }     
        }
    }
    double dt2;
    MPI_Allreduce( &dt_m , &dt2 , 1 , MPI_DOUBLE , MPI_MIN , sim_comm );

    //double dt_newt = cell_mindt_newt(theCells, theSim);
    //printf("\nNewtonian dt: %lg, GR dt: %lg\n", dt_newt, dt2);

    return dt2;
}

//TODO: WRITE THIS
double cell_mindt_grdisc(struct Cell ***theCells, struct Sim *theSim)
{
    return 0.0;
}
