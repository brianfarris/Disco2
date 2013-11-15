#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Sim.h"
#include "../Headers/Face.h"
#include "../Headers/GravMass.h"
#include "../Headers/Metric.h"
#include "../Headers/header.h"

double fgrav( double M , double r , double eps, double n ){
  return( M*pow(r,n-1.)/pow( pow(r,n) + pow(eps,n) ,1.+1./n) );
}

double fgrav_neg_centrifugal( double M , double r , double eps, double n ){
  double Om = 20.0;
  return( M*r*Om*Om );
}

void gravMassForce( struct GravMass * theGravMasses ,struct Sim * theSim, int p , double r , double phi , double * fr , double * fp ){
  double G_EPS=sim_G_EPS(theSim);
  double PHI_ORDER=sim_PHI_ORDER(theSim);

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

void cell_add_src( struct Cell *** theCells ,struct Sim * theSim, struct GravMass * theGravMasses , double dt ){
  int GRAV2D=sim_GRAV2D(theSim);
  int i,j,k;
  //FILE * gradpsifile= fopen("gradpsi.dat","w");
  for( k=0 ; k<sim_N(theSim,Z_DIR) ; ++k ){
    double zm = sim_FacePos(theSim,k-1,Z_DIR);
    double zp = sim_FacePos(theSim,k,Z_DIR);
    for( i=0 ; i<sim_N(theSim,R_DIR) ; ++i ){
      double rm = sim_FacePos(theSim,i-1,R_DIR);
      double rp = sim_FacePos(theSim,i,R_DIR);
      for( j=0 ; j<sim_N_p(theSim,i) ; ++j ){
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
        for( p=0 ; p<sim_NumGravMass(theSim); ++p ){
          gravMassForce( theGravMasses , theSim , p , gravdist , phi , &fr , &fp );
          Fr += fr;
          Fp += fp;
        }

        if(sim_Background(theSim) == NEWTON)
        {
            c->cons[SRR] += dt*dV*( rho*vp*vp + Pp )/r;
            c->cons[SRR] += dt*dV*rho*Fr*sint;
            c->cons[LLL] += dt*dV*rho*Fp*r;
            c->cons[SZZ] += dt*dV*rho*Fr*cost;
            c->cons[TAU] += dt*dV*rho*( Fr*(vr*sint+vz*cost) + Fp*vp );
        }
        else if(sim_Background(theSim) == GR)
        {
            int i,j,k;
            double a, b[3], sqrtg, n[4], u[4], u_d[4], s, sk, rhoh, GAMMALAW;
            double v[3];
            struct Metric *g;
            
            g = metric_create(time_global, r, phi, z);
            a = metric_lapse(g);
            for(i=0; i<3; i++)
                b[i] = metric_shift_u(g,i);
            sqrtg = metric_sqrtgamma(g)/r;

            v[0] = c->prim[URR];
            v[1] = c->prim[UPP];
            v[2] = c->prim[UZZ];
            //Contravariant Four-Velocity u[i] = u^i
            u[0] = 1.0 / sqrt(-metric_g_dd(g,0,0) - 2*metric_dot3_u(g,b,v) - metric_square3_u(g,v));
            u[1] = u[0] * c->prim[URR];
            u[2] = u[0] * c->prim[UPP];
            u[3] = u[0] * c->prim[UZZ];
            //Covariant Four-Velocity u_d[i] = u_i
            for(i=0; i<4; i++)
                for(j=0; j<4; j++)
                    u_d[i] = metric_g_uu(g,i,j) * u[j];
            //Normal vector n[i] = n^i
            n[0] = 1.0/a;
            n[1] = -n[0]*b[0];
            n[2] = -n[0]*b[1];
            n[3] = -n[0]*b[2];

            GAMMALAW = sim_GAMMALAW(theSim);
            rhoh = rho + GAMMALAW*Pp/(GAMMALAW-1);
            

            //Momentum sources and contribution to energy source
            s = 0;
            for(k=0; k<4; k++)
                if(!metric_killcoord(g,k))
                {
                    sk = 0;
                    for(i=0; i<4; i++)
                    {
                        sk += 0.5*(rhoh*u[i]*u[i]+Pp*metric_g_uu(g,i,i)) * metric_dg_dd(g,k,i,i);
                        for(j=i+1; j<4; j++)
                            sk += (rhoh*u[i]*u[j]+Pp*metric_g_uu(g,i,j)) * metric_dg_dd(g,k,i,j);
                    }
                    if(k == 1)
                    {
                        c->cons[SRR] += dt*dV*sqrtg*a * sk;
                        /*
                        printf("sqrtg = %lg, a = %lg, rhoh = %lg\n", sqrtg, a, rhoh);
                        printf("u0 = %lg, ur = %lg, up = %lg, uz = %lg\n", u[0], u[1], u[2], u[3]);
                        printf("    (%lg %lg %lg %lg)\n", metric_g_uu(g,0,0), metric_g_uu(g,0,1), metric_g_uu(g,0,2), metric_g_uu(g,0,3));
                        printf("g = (%lg %lg %lg %lg)\n", metric_g_uu(g,1,0), metric_g_uu(g,1,1), metric_g_uu(g,1,2), metric_g_uu(g,1,3));
                        printf("    (%lg %lg %lg %lg)\n", metric_g_uu(g,2,0), metric_g_uu(g,2,1), metric_g_uu(g,2,2), metric_g_uu(g,2,3));
                        printf("    (%lg %lg %lg %lg)\n", metric_g_uu(g,3,0), metric_g_uu(g,3,1), metric_g_uu(g,3,2), metric_g_uu(g,3,3));
                        printf("     (%lg %lg %lg %lg)\n", metric_dg_dd(g,1,0,0), metric_dg_dd(g,1,0,1), metric_dg_dd(g,1,0,2), metric_dg_dd(g,1,0,3));
                        printf("dg = (%lg %lg %lg %lg)\n", metric_dg_dd(g,1,1,0), metric_dg_dd(g,1,1,1), metric_dg_dd(g,1,1,2), metric_dg_dd(g,1,1,3));
                        printf("     (%lg %lg %lg %lg)\n", metric_dg_dd(g,1,2,0), metric_dg_dd(g,1,2,1), metric_dg_dd(g,1,2,2), metric_dg_dd(g,1,2,3));
                        printf("     (%lg %lg %lg %lg)\n", metric_dg_dd(g,1,3,0), metric_dg_dd(g,1,3,1), metric_dg_dd(g,1,3,2), metric_dg_dd(g,1,3,3));
                        printf("Source Sr (r=%lg, rho=%lg, Pp=%lg, vr=%lg, vp=%lg): %lg\n", r, rho, Pp, v[0], v[1], sqrtg*a*sk);
                        */
                    }
                    else if(k == 2)
                    {
                        c->cons[LLL] += dt*dV*sqrtg*a * sk;
                        printf("Source Sp: %lg\n", sqrtg*a*sk);
                    }
                    else if(k == 3)
                    {
                        c->cons[SZZ] += dt*dV*sqrtg*a * sk;
                        printf("Source Sz: %lg\n", sqrtg*a*sk);
                    }
                    s -= n[k]*sk;
                }
            //Remaining energy sources
            for(k=0; k<4; k++)
                if(!metric_killcoord(g,k))
                {
                    sk = (rhoh*u[0]*u[k] + Pp*metric_g_uu(g,0,k))*metric_dlapse(g,k);
                    for(i=0; i<4; i++)
                    {
                        if(i == k)
                            sk += a*(rhoh*u[0]*u_d[i]+Pp)*metric_dg_uu(g,k,0,i);
                        else
                            sk += a*rhoh*u[0]*u_d[i]*metric_dg_uu(g,k,0,i);

                    }
                }
            s += sk;

            c->cons[TAU] += dt*dV*sqrtg*a * s;

            metric_destroy(g);
        }
        else
        {
            printf("ERROR: Unknown Background in cell_add_src().\n");
            exit(0);
        }
      }
    }
  }
}

