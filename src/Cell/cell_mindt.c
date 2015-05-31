#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/EOS.h"
#include "../Headers/Metric.h"
#include "../Headers/Sim.h"
#include "../Headers/header.h"

//Local functions
double cell_mindt_binary_step(double r, double phi, double d4V, double tau,
                                double rhoh, double *u, double *U,
                                struct Sim *theSim);

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
        vm = fabs((vn - sig*bn - dv) / (1.0+sig) - w/r);
        vp = fabs((vn - sig*bn + dv) / (1.0+sig) - w/r);
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

double cell_maxvel_grdisc(double *prim, int dir, double w, double r, double phi, double z, struct Sim *theSim)
{
    double maxv;
    double rho, Pp, v[3], eps, rhoh, vn;
    double a, bn, b[3], u02;
    double vp, vm, sig, dv, gam, cs2;
    int i;
    struct Metric *g;

    rho = prim[RHO];
    v[0] = prim[URR];
    v[1] = prim[UPP];
    v[2] = prim[UZZ];
    vn = v[dir];

    Pp = eos_ppp(prim, theSim);
    eps = eos_eps(prim, theSim);
    cs2 = eos_cs2(prim, theSim);
    rhoh = rho + rho*eps + Pp;
    
    g = metric_create(time_global, r, phi, z, theSim);
    a = metric_lapse(g);
    for(i=0; i<3; i++)
        b[i] = metric_shift_u(g,i);
    bn = b[dir];
    gam = metric_gamma_uu(g,dir,dir);
    u02 = 1.0/(-metric_g_dd(g,0,0)-2.0*metric_dot3_u(g,b,v)-metric_square3_u(g,v));

    sig = cs2 / ((1-cs2)*u02*a*a);
    dv = sqrt(sig*(1.0+sig)*a*a*gam - sig*(vn+bn)*(vn+bn));

    if(dir == PDIRECTION)
    {
        vm = fabs((vn - sig*bn - dv) / (1.0+sig) - w/r);
        vp = fabs((vn - sig*bn + dv) / (1.0+sig) - w/r);
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
            
            //if(r < metric_horizon(theSim))
            //    continue;

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
                    int mu;
                    double rho = theCells[k][i][j].prim[RHO];
                    double Pp = theCells[k][i][j].prim[PPP];
                    double GAM = sim_GAMMALAW(theSim);
                    double u0, b[3], v[3];
                    double rhoh = rho + GAM/(GAM-1)*Pp;
                    double nu;
                    double alpha = sim_AlphaVisc(theSim);
                    struct Metric *g;

                    v[0] = theCells[k][i][j].prim[URR];
                    v[1] = theCells[k][i][j].prim[UPP];
                    v[2] = theCells[k][i][j].prim[UZZ];
                    g = metric_create(time_global, r, phi, z, theSim);
                    for(mu=0; mu<3; mu++)
                        b[mu] = metric_shift_d(g,mu);
                    u0 = 1.0 / sqrt(-metric_g_dd(g,0,0)-metric_square3_u(g,v)
                                    - 2*(v[0]*b[0]+v[1]*b[1]+v[2]*b[2]));
                    metric_destroy(g);
                    double H = sqrt(r*r*r*Pp / (M*rhoh)) / u0;

                    if(alpha < 0.0)
                        nu = -alpha;
                    else
                    {
                        nu = alpha * sqrt(GAM*Pp/rhoh) * H;
                    }
                    double dtv = 0.25*dx*dx/nu;
                    if(alpha == 0.0 || rp < metric_horizon(theSim))
                        dtv = 1.0e100; //HUGE
                    
                    //Luminosity Time
                    double maxdisp = 0.1; //Should be a fraction like 0.1
                    double dtl;

                    //TODO: Decide whether this should be here at all.  If not, remove it.
                    if(0 && sim_CoolingType(theSim) != COOL_NONE)
                    {
                        double dV = 0.5*(rp+rm)*dr*dphi*dz;
                        double U[4];
                        double cool = 0.0;

                        if(rp > metric_horizon(theSim))
                            cool = eos_cool(theCells[k][i][j].prim, H, theSim);

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
                        if(dtl < dtl_min && dtl > 0.0)
                            dtl_min = dtl;
                        dtl = maxdisp * fabs(theCells[k][i][j].cons[LLL] / (dV*al*sqrtg*u_d[2]*cool));
                        if(dtl < dtl_min && dtl > 0.0)
                            dtl_min = dtl;
                        if(sim_N(theSim, Z_DIR) > 1)
                        {
                            dtl = maxdisp * fabs(theCells[k][i][j].cons[SZZ] / (dV*al*sqrtg*u_d[3]*cool));
                            if(dtl < dtl_min && dtl > 0.0)
                                dtl_min = dtl;
                        }
                        dtl = dtl_min;

                        if(dtl == 0.0)
                            printf("r=%.12g (%.12g %.12g %.12g)\n",r, dt, dtv, dtl);

                        dt = dt/(1.0 + dt/dtv + dt/dtl);
                    }
                    else
                    {
                        dt = dt/(1.0 + dt/dtv);
                    }

                }

                if(sim_BoostType(theSim) == BOOST_BIN)
                {
                    struct Metric *g;
                    g = metric_create(time_global, r, phi, z, theSim);
                    double al, b[3], sqrtg, v[3], u[4], U[4];
                    al = metric_lapse(g);
                    b[0] = metric_shift_u(g,0);
                    b[1] = metric_shift_u(g,1);
                    b[2] = metric_shift_u(g,2);
                    sqrtg = metric_sqrtgamma(g)/r;
                    v[0] = theCells[k][i][j].prim[URR];
                    v[1] = theCells[k][i][j].prim[UPP];
                    v[2] = theCells[k][i][j].prim[UZZ];
                    u[0] = 1.0/sqrt(-metric_g_dd(g,0,0) - 2*metric_dot3_u(g,b,v) - metric_square3_u(g,v));
                    u[1] = u[0]*v[0]; u[2] = u[0]*v[1]; u[3] = u[0]*v[2];
                    U[0] = metric_frame_U_u(g, 0, theSim);
                    U[1] = metric_frame_U_u(g, 1, theSim);
                    U[2] = metric_frame_U_u(g, 2, theSim);
                    U[3] = metric_frame_U_u(g, 3, theSim);
                    metric_destroy(g);
                    double rho = theCells[k][i][j].prim[RHO];
                    double Pp = theCells[k][i][j].prim[PPP];
                    double GAM = sim_GAMMALAW(theSim);
                    double rhoh = rho + GAM/(GAM-1)*Pp;
                    double dV = 0.5*(rp+rm)*dr*dphi*dz;
                    double dtb = cell_mindt_binary_step(r, phi, al*sqrtg*dV, 
                                    theCells[k][i][j].cons[TAU], rhoh, u, U,
                                    theSim);
                    dt = dt/(1.0 + dt/dtb);
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

double cell_mindt_grdisc(struct Cell ***theCells, struct Sim *theSim)
{
    double dt_m = 1.e100;//HUGE_VAL;

    double M = sim_GravM(theSim);
    double alpha = sim_AlphaVisc(theSim);
    double rho, eps, Pp, rhoh, cs2, H, nu_visc, qdot;
    double al, b[3], sqrtg, v[3], u[4], u_d[4], U[4];
    struct Metric *g;
    int i,j,k, mu;

    for(k = sim_Nghost_min(theSim,Z_DIR); k < sim_N(theSim,Z_DIR)-sim_Nghost_max(theSim,Z_DIR); k++)
    {
        double zm = sim_FacePos(theSim,k-1,Z_DIR);
        double zp = sim_FacePos(theSim,k,Z_DIR);
        double dz = zp-zm;
        double z = 0.5*(zm+zp);
        
        for(i = sim_Nghost_min(theSim,R_DIR); i < sim_N(theSim,R_DIR)-sim_Nghost_max(theSim,R_DIR); i++)
        {
            double rm = sim_FacePos(theSim,i-1,R_DIR);
            double rp = sim_FacePos(theSim,i,R_DIR);
            double dr = rp-rm;
            double r = .5*(rp+rm);
            
            //TODO: cleaner way of doing this?
            if(r < metric_horizon(theSim))
                continue;

            for(j = 0; j < sim_N_p(theSim,i); j++)
            {
                int jm = j-1;
                if(j == 0) 
                    jm = sim_N_p(theSim,i)-1;
                double w = .5*(theCells[k][i][j].wiph+theCells[k][i][jm].wiph); 
                double dphi = theCells[k][i][j].dphi;
                double phi = theCells[k][i][j].tiph - 0.5*dphi;
                double dV = 0.5*(rp+rm)*dr*dphi*dz;

                double ar, ap, az;
                ar = cell_maxvel_grdisc(theCells[k][i][j].prim, 0, w, r, phi, 
                                        z, theSim);
                ap = cell_maxvel_grdisc(theCells[k][i][j].prim, 1, w, r, phi, 
                                        z, theSim);
                az = cell_maxvel_grdisc(theCells[k][i][j].prim, 2, w, r, phi, 
                                        z, theSim);
                
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

                v[0] = theCells[k][i][j].prim[URR];
                v[1] = theCells[k][i][j].prim[UPP];
                v[2] = theCells[k][i][j].prim[UZZ];

                g = metric_create(time_global, r, phi, z, theSim);

                for(mu=0; mu<3; mu++)
                    b[mu] = metric_shift_d(g,mu);
                u[0] = 1.0 / sqrt(-metric_g_dd(g,0,0) - metric_square3_u(g,v)
                                - 2*(v[0]*b[0]+v[1]*b[1]+v[2]*b[2]));
                u[1] = u[0]*v[0];
                u[2] = u[0]*v[1];
                u[3] = u[0]*v[2];

                rho = theCells[k][i][j].prim[RHO];
                Pp = eos_ppp(theCells[k][i][j].prim, theSim);
                eps = eos_eps(theCells[k][i][j].prim, theSim);
                cs2 = eos_cs2(theCells[k][i][j].prim, theSim);

                rhoh = rho + rho*eps + Pp;
                if(sim_Metric(theSim) == KERR_KS)
                {
                    double A = M*sim_GravA(theSim);
                    double l2 = 0.0;
                    double e2 = 0.0;
                    int mu;

                    for(mu = 0; mu < 4; mu++)
                    {
                        e2 += metric_g_dd(g,0,mu)*u[mu];
                        l2 += metric_g_dd(g,2,mu)*u[mu];
                    }
                    e2 = e2*e2;
                    l2 = l2*l2;

                    H = r*r * sqrt(Pp / (rhoh * (l2-A*A*(e2-1.0))));

                }
                else
                    H = sqrt(r*r*r*Pp / (M*rhoh)) / u[0];

                if(alpha < 0.0)
                    nu_visc = -alpha;
                else
                    nu_visc = alpha * sqrt(cs2) * H;
                double dtv = 0.25*dx*dx/nu_visc;
                if(alpha == 0.0)
                    dtv = 1.0e100; //HUGE
                
                //Luminosity Time
                double maxdisp = 0.1; //Should be a fraction like 0.1
                double dtl;
                
                //TODO: Remove?
                if(0 && (sim_CoolingType(theSim) != COOL_NONE))
                {
                    qdot = eos_cool(theCells[k][i][j].prim, H, theSim);
                    
                    al = metric_lapse(g);
                    sqrtg = metric_sqrtgamma(g)/r;
                    for(mu = 0; mu < 4; mu++)
                        U[mu] = metric_frame_U_u(g, mu, theSim);
                    
                    for(mu=0; mu<4; mu++)
                        u_d[mu] = 0.0;
                    for(mu=0; mu<4; mu++)
                    {
                        u_d[0] += metric_g_dd(g,0,mu)*u[mu];
                        u_d[1] += metric_g_dd(g,1,mu)*u[mu];
                        u_d[2] += metric_g_dd(g,2,mu)*u[mu];
                        u_d[3] += metric_g_dd(g,3,mu)*u[mu];
                    }

                    double dtl_min;
                    dtl = maxdisp * fabs(theCells[k][i][j].cons[TAU]
                                        / (dV*sqrtg*al* qdot
                        * (u_d[0]*U[0]+u_d[1]*U[1]+u_d[2]*U[2]+u_d[3]*U[3])));
                    dtl = fabs(dtl);
                    if (dtl < 0)
                        printf("WHAT.");
                    dtl_min = dtl;
                    dtl = maxdisp * fabs(theCells[k][i][j].cons[SRR]
                                        / (dV*al*sqrtg*u_d[1]*qdot));
                    if(dtl < dtl_min)
                        dtl_min = dtl;
                    dtl = maxdisp * fabs(theCells[k][i][j].cons[LLL]
                                        / (dV*al*sqrtg*u_d[2]*qdot));
                    if(dtl < dtl_min)
                        dtl_min = dtl;
                    if(sim_N(theSim, Z_DIR) > 1)
                    {
                        dtl = maxdisp * fabs(theCells[k][i][j].cons[SZZ]
                                            / (dV*al*sqrtg*u_d[3]*qdot));
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

                if(dt_m > dt)
                    dt_m = dt;

                metric_destroy(g);
            }     
        }
    }
    double dt2;
    MPI_Allreduce( &dt_m , &dt2 , 1 , MPI_DOUBLE , MPI_MIN , sim_comm );

    //double dt_newt = cell_mindt_newt(theCells, theSim);
    //printf("\nNewtonian dt: %lg, GR dt: %lg\n", dt_newt, dt2);

    return dt2;
}

double cell_mindt_binary_step(double r, double phi, double d4V, double tau,
                                double rhoh, double *u, double *U, 
                                struct Sim *theSim)
{
    //Calculates the time-scale for energy loss due to the binary.
    
    double a = sim_BinA(theSim);
    double M2 = sim_BinM(theSim);

    double Fr, Fp, Ft;
    double X = 2*a + r*cos(phi);
    double Y = r*sin(phi);
    double R = sqrt(X*X+Y*Y);

    Fr = -M2/(R*R)*( X*cos(phi)/R + Y*sin(phi)/R)*rhoh*u[0]*u[0];
    Fp = -M2/(R*R)*(-X*sin(phi)/R + Y*cos(phi)/R)*rhoh*u[0]*u[0]*r;
    Ft = -(u[1]*Fr + u[2]*Fp/r)/u[0];

    return 0.5 * fabs(tau / (d4V * (-U[0]*Ft-U[1]*Fr-U[2]*Fp)));
}
