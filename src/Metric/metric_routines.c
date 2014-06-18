#define METRIC_PRIVATE_DEFS
#include <stdio.h>
#include <math.h>
#include "../Headers/Metric.h"
#include "../Headers/Sim.h"
#include "../Headers/header.h"


double metric_square3_u(struct Metric *g, double *vec)
{
    int i,j;
    double v2 = 0;
    for(i=0; i<3; i++)
        for(j=0; j<3; j++)
            v2 += metric_gamma_dd(g,i,j)*vec[i]*vec[j];
    return v2;
}

double metric_square3_d(struct Metric *g, double *vec)
{
    int i,j;
    double v2 = 0.0;
    for(i=0; i<3; i++)
        for(j=0; j<3; j++)
            v2 += metric_gamma_uu(g,i,j)*vec[i]*vec[j];
    return v2;
}

double metric_square4_u(struct Metric *g, double *vec)
{
    int i,j;
    double v2 = 0;
    for(i=0; i<4; i++)
        for(j=0; j<4; j++)
            v2 += metric_g_dd(g,i,j)*vec[i]*vec[j];
    return v2;
}

double metric_square4_d(struct Metric *g, double *vec)
{
    int i,j;
    double v2 = 0;
    for(i=0; i<4; i++)
        for(j=0; j<4; j++)
            v2 += metric_g_uu(g,i,j)*vec[i]*vec[j];
    return v2;
}

double metric_dot3_u(struct Metric *g, double *a, double *b)
{
    int i,j;
    double v2 = 0;
    for(i=0; i<3; i++)
        for(j=0; j<3; j++)
            v2 += metric_gamma_dd(g,i,j)*a[i]*b[j];
    return v2;
}

double metric_dot3_d(struct Metric *g, double *a, double *b)
{
    int i,j;
    double v2 = 0;
    for(i=0; i<3; i++)
        for(j=0; j<3; j++)
            v2 += metric_gamma_uu(g,i,j)*a[i]*b[j];
    return v2;
}

double metric_dot4_u(struct Metric *g, double *a, double *b)
{
    int i,j;
    double v2 = 0;
    for(i=0; i<4; i++)
        for(j=0; j<4; j++)
            v2 += metric_g_dd(g,i,j)*a[i]*b[j];
    return v2;
}

double metric_dot4_d(struct Metric *g, double *a, double *b)
{
    int i,j;
    double v2 = 0;
    for(i=0; i<4; i++)
        for(j=0; j<4; j++)
            v2 += metric_g_uu(g,i,j)*a[i]*b[j];
    return v2;
}

double metric_conn(struct Metric *g, int tau, int mu, int nu)
{
    //Returns value of the connection coefficient \Gamma^\tau_{\mu\nu}
    double connection = 0.0;
    int sig;
    for(sig=0; sig<4; sig++)
        connection += metric_g_uu(g,tau,sig)*(metric_dg_dd(g,mu,sig,nu) + metric_dg_dd(g,nu,mu,sig) - metric_dg_dd(g,sig,mu,nu));
    return 0.5*connection;
}

void metric_shear_uu(struct Metric *g, double *v, double *dv, double *shear, struct Sim *theSim)
{
    //Fills shear tensor sigma^{ij} = shear(4*i+j) given 
    // coordinate 3-velocity v and (partial der) gradient dv (d_{i}v^{j} =  dv[3*i+j])

    double divu = 0.0;
    double cdu, val, du0;
    double a, b[3], u[4], du[16];
    int mu,nu,la;

    for(mu=0; mu<16; mu++)
    {
        du[mu] = 0.0;
        shear[mu] = 0.0;
    }

    a = metric_lapse(g);
    for(mu=0; mu<3; mu++)
        b[mu] = metric_shift_u(g, mu);

    u[0] = 1.0/sqrt(-metric_g_dd(g,0,0) - 2*metric_dot3_u(g,b,v)-metric_square3_u(g,v));
    u[1] = u[0]*v[0];
    u[2] = u[0]*v[1];
    u[3] = u[0]*v[2];

    for(mu=0; mu<4; mu++)
    {
        du0 = metric_dg_dd(g,mu,0,0);
        for(nu=1; nu<4; nu++)
        {
            du0 += 2*v[nu-1]*metric_dg_dd(g,mu,0,nu) + 2*metric_g_dd(g,0,nu)*dv[3*mu+nu-1];
            for(la=1; la<4; la++)
                du0 += v[nu-1]*v[la-1]*metric_dg_dd(g,mu,nu,la) + 2*metric_g_dd(g,nu,la)*v[nu-1]*dv[3*mu+nu-1];
        }
        du[4*mu] = 0.5*u[0]*u[0]*u[0]*du0;
    }
    
    for(mu=1; mu<4; mu++)
        for(nu=1; nu<4; nu++)
            du[4*mu+nu] = v[nu-1]*du[4*mu] + u[0]*dv[3*mu+nu-1];

    for(mu=0; mu<4; mu++)
    {
        for(nu=0; nu<4; nu++)
        {
            cdu = du[4*mu+nu];  //d_mu u^nu
            for(la=0; la<4; la++)
                cdu += metric_conn(g, nu,mu,la)*u[la];

            for(la=0; la<4; la++)
            {
                val = cdu*(metric_g_uu(g,mu,la)+u[mu]*u[la]);
                shear[4*nu+la] += val;
                shear[4*la+nu] += val;
            }
            if(mu == nu)
                divu += cdu;
        }
    }

    if(sim_N_global(theSim, Z_DIR) != 1)
        divu *= 2.0/3.0;

    for(mu=0; mu<4; mu++)
        for(nu=0; nu<4; nu++)
            shear[4*mu+nu] -= divu*(metric_g_uu(g,mu,nu)+u[mu]*u[nu]);
}

