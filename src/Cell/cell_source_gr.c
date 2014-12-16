#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/EOS.h"
#include "../Headers/Sim.h"
#include "../Headers/Face.h"
#include "../Headers/GravMass.h"
#include "../Headers/Metric.h"
#include "../Headers/header.h"

void cell_add_src_gr( struct Cell *** theCells ,struct Sim * theSim, struct GravMass * theGravMasses , double dt )
{
    int i,j,k;
  
    FILE *sourcefile;
    if(PRINTTOOMUCH)
        sourcefile = fopen("source.out","w");
  
    for( k=0 ; k<sim_N(theSim,Z_DIR) ; ++k )
    {
        double zm = sim_FacePos(theSim,k-1,Z_DIR);
        double zp = sim_FacePos(theSim,k,Z_DIR);
    
        for( i=0 ; i<sim_N(theSim,R_DIR) ; ++i )
        {
            double rm = sim_FacePos(theSim,i-1,R_DIR);
            double rp = sim_FacePos(theSim,i,R_DIR);

            for( j=0 ; j<sim_N_p(theSim,i) ; ++j )
            {
                struct Cell * c = &(theCells[k][i][j]);
                double phi = c->tiph-.5*c->dphi;
                double dphi = c->dphi;

                double rho = c->prim[RHO];
                double Pp  = c->prim[PPP];
                double r   = .5*(rp+rm);
                double z  = .5*(zp+zm);
                double dz = zp-zm;
                double dV = dphi*.5*(rp*rp-rm*rm)*dz;

       
                int mu, nu, la;
                double a, b[3], sqrtg, u[4], u_d[4], U[4], s, sk, rhoh, GAMMALAW;
                double v[3];
                struct Metric *g;
                
                g = metric_create(time_global, r, phi, z, theSim);
                a = metric_lapse(g);
                for(mu=0; mu<3; mu++)
                    b[mu] = metric_shift_u(g,mu);
                sqrtg = metric_sqrtgamma(g)/r;
                for(mu=0; mu<4; mu++)
                    U[mu] = metric_frame_U_u(g, mu, theSim);

                v[0] = c->prim[URR];
                v[1] = c->prim[UPP];
                v[2] = c->prim[UZZ];
                //Contravariant Four-Velocity u[i] = u^i
                u[0] = 1.0 / sqrt(-metric_g_dd(g,0,0) - 2*metric_dot3_u(g,b,v) - metric_square3_u(g,v));
                u[1] = u[0] * c->prim[URR];
                u[2] = u[0] * c->prim[UPP];
                u[3] = u[0] * c->prim[UZZ];
                //Covariant Four-Velocity u_d[i] = u_i
                for(mu=0; mu<4; mu++)
                {
                    u_d[mu] = 0.0;
                    for(nu=0; nu<4; nu++)
                        u_d[mu] += metric_g_dd(g,mu,nu) * u[nu];
                }

                GAMMALAW = sim_GAMMALAW(theSim);
                rhoh = rho + GAMMALAW*Pp/(GAMMALAW-1);
               
                //Viscous Terms
                double cool, visc[16], viscm[16];
                double M = sim_GravM(theSim);
                cool = 0.0;
                for(mu=0; mu<16; mu++)
                {
                    visc[mu] = 0.0;
                    viscm[mu] = 0.0;
                }
                if(sim_Background(theSim) == GRVISC1)
                {
                    double dv[12];
                    double cs2;
                    double alpha = sim_AlphaVisc(theSim);

                    double height = sqrt(Pp*r*r*r/(rhoh*M))/u[0];

                    cs2 = GAMMALAW*Pp / rhoh;
                    if (alpha > 0)
                        alpha *= sqrt(cs2) * height * rho;
                    else
                        alpha = -alpha * rho;

                    dv[0] = 0.0;  dv[1] = 0.0;  dv[2] = 0.0;
                    dv[3] = cell_gradr(c, URR);  dv[4] = cell_gradr(c, UPP);  dv[5] = cell_gradr(c, UZZ);
                    dv[6] = cell_gradp(c, URR);  dv[7] = cell_gradp(c, UPP);  dv[8] = cell_gradp(c, UZZ);
                    dv[9] = cell_gradz(c, URR);  dv[10] = cell_gradz(c, UPP);  dv[11] = cell_gradz(c, UZZ);
                    
                    metric_shear_uu(g, v, dv, visc, theSim);
                    if(sim_N(theSim,Z_DIR)==1)
                    {
                        for(mu=0; mu<4; mu++)
                        {
                            visc[4*mu+3] = 0.0;
                            visc[4*3+mu] = 0.0;
                        }
                    }
                
                    for(mu=0; mu<16; mu++)
                        visc[mu] *= -alpha;
                    //viscm[4*mu+nu] = visc^mu_nu
                    for(mu=0; mu<4; mu++)
                        for(nu=0; nu<4; nu++)
                            for(la=0; la<4; la++)
                                viscm[4*mu+nu] += metric_g_dd(g,nu,la)*visc[4*mu+la];
                  
                    cool = eos_cool(c->prim, height, theSim);
                }

                //Momentum sources and contribution to energy source
                s = 0;
                int max_dim = 4;
                if(sim_N(theSim,Z_DIR)==1)
                    max_dim = 3;
       
                for(la=0; la<max_dim; la++)
                    if(!metric_killcoord(g,la))
                    {
                        sk = 0;
                        for(mu=0; mu<max_dim; mu++)
                        {
                            sk += 0.5*(rhoh*u[mu]*u[mu]+Pp*metric_g_uu(g,mu,mu)+visc[4*mu+mu]) * metric_dg_dd(g,la,mu,mu);
                            for(nu=mu+1; nu<max_dim; nu++)
                                sk += (rhoh*u[mu]*u[nu]+Pp*metric_g_uu(g,mu,nu)+visc[4*mu+nu]) * metric_dg_dd(g,la,mu,nu);
                        }
                        if(la == 1)
                        {
                            c->cons[SRR] += dt*dV*sqrtg*a * sk;
                            if(PRINTTOOMUCH)
                            {
                                printf("SRR source: (%d,%d,%d): r=%.12g, dV=%.12g, s = %.12g, S = %.12g\n",i,j,k,r,dV,a*sqrtg*sk,dt*dV*sqrtg*a * sk);
                                //fprintf(sourcefile, "(%d,%d,%d): r=%.12g, Fg=%.12g, Fc=%.12g, P/r=%.12g, F=%.12g\n", 
                                //        i,j,k,r,0.5*rhoh*u[0]*u[0]*metric_dg_dd(g,1,0,0), 0.5*rhoh*u[2]*u[2]*metric_dg_dd(g,1,2,2), 
                                //        0.5*metric_g_uu(g,2,2)*Pp*metric_dg_dd(g,1,2,2), sk);
                            }
                        }
                        else if(la == 2)
                        {
                            c->cons[LLL] += dt*dV*sqrtg*a * sk;
                        }
                        else if(la == 3)
                        {
                            c->cons[SZZ] += dt*dV*sqrtg*a * sk;
                            if(PRINTTOOMUCH)
                            {
                                printf("SZZ source: (%d,%d,%d): r=%.12f, dV=%.12f, s = %.12f, S = %.12f\n",i,j,k,r,dV,a*sqrtg*sk,dt*dV*sqrtg*a * sk);
                            }
                        }
                        s -= U[la]*sk;
                    }

                //Remaining energy sources
                
                for(mu=0; mu<max_dim; mu++)
                    for(nu=0; nu<max_dim; nu++)
                    {
                        if(mu == nu)
                            s -= (rhoh*u[mu]*u_d[nu] + Pp + viscm[4*mu+nu]) * metric_frame_dU_du(g,mu,nu,theSim);
                        else
                            s -= (rhoh*u[mu]*u_d[nu] + viscm[4*mu+nu]) * metric_frame_dU_du(g,mu,nu,theSim);
                    }
                
                if(PRINTTOOMUCH)
                {
                    printf("TAU source: (%d,%d,%d): r=%.12g, dV=%.12g, s = %.12g, S = %.12g\n",i,j,k,r,dV,a*sqrtg*s,dt*dV*sqrtg*a * s);
                }

                c->cons[TAU] += dt*dV*sqrtg*a * s;

                //Cooling
                //Cooling should never (?) change the sign of the momenta,
                //if it will, instead reduce the momenta by 0.9.  These lines
                //should be irrelevant because of the luminosity-limited timestep.
                if(fabs(c->cons[SRR]) < fabs(dt*dV*sqrtg*a* cool*u_d[1]))
                    c->cons[SRR] /= 1.1;
                else
                    c->cons[SRR] -= dt*dV*sqrtg*a* cool*u_d[1];
                if(fabs(c->cons[LLL]) <  fabs(dt*dV*sqrtg*a* cool*u_d[2]))
                    c->cons[LLL] /= 1.1;
                else
                    c->cons[LLL] -= dt*dV*sqrtg*a* cool*u_d[2];
                if( fabs(c->cons[SZZ]) < fabs(dt*dV*sqrtg*a* cool*u_d[3]))
                    c->cons[SZZ] /= 1.1;
                else
                    c->cons[SZZ] -= dt*dV*sqrtg*a* cool*u_d[3];
                if(fabs(c->cons[TAU]) < fabs(dt*dV*sqrtg*a* cool * (u_d[0]*U[0]+u_d[1]*U[1]+u_d[2]*U[2]+u_d[3]*U[3])))
                {
                    c->cons[TAU] /= 1.1;
                    printf("That's pretty cool!\n");
                }
                else
                    c->cons[TAU] += dt*dV*sqrtg*a* cool * (u_d[0]*U[0]+u_d[1]*U[1]+u_d[2]*U[2]+u_d[3]*U[3]); 
                if(PRINTTOOMUCH)
                {
                    printf("SRR cooling: (%d,%d,%d): r=%.12g, dV=%.12g, q = %.12g, Q = %.12g\n",i,j,k,r,dV,-a*sqrtg* cool*u_d[1],-dt*dV*sqrtg*a *cool*u_d[1]);
                    printf("LLL cooling: (%d,%d,%d): r=%.12g, dV=%.12g, q = %.12g, Q = %.12g\n",i,j,k,r,dV,-a*sqrtg* cool*u_d[2],-dt*dV*sqrtg*a *cool*u_d[2]);
                    printf("TAU cooling: (%d,%d,%d): r=%.12g, dV=%.12g, q = %.12g, Q = %.12g\n",i,j,k,r,dV,a*sqrtg* cool * (u_d[0]*U[0]+u_d[1]*U[1]+u_d[2]*U[2]+u_d[3]*U[3]),dt*dV*sqrtg*a *cool * (u_d[0]*U[0]+u_d[1]*U[1]+u_d[2]*U[2]+u_d[3]*U[3]));
                }

                metric_destroy(g);
            }
        }
    }

    if(PRINTTOOMUCH)
        fclose(sourcefile);
}

