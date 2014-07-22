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

double metric_frame_U_u_euler(struct Metric *g, int mu, struct Sim *theSim)
{
    if(mu == 0)
        return 1.0/metric_lapse(g);
    else if (mu > 0 && mu < 4)
        return metric_shift_u(g,mu-1) / metric_lapse(g);
    return 0.0;
}

double metric_frame_dU_du_euler(struct Metric *g, int mu, int nu, struct Sim *theSim)
{
    if(nu == 0)
    {
        double a = metric_lapse(g);
        return -metric_dlapse(g,mu) / (a*a);
    }
    else if (nu > 0 && nu < 4)
    {
        double a = metric_lapse(g);
        return (metric_dlapse(g,mu)*metric_shift_u(g,nu-1) - a*metric_dshift_u(g,mu,nu-1)) / (a*a);
    }
    return 0.0;
}

double metric_frame_U_u_kep(struct Metric *g, int mu, struct Sim *theSim)
{
    double M = sim_GravM(theSim);
    double r = g->x[1];

    double u0, vr, vp;

    if(sim_Metric(theSim) == SCHWARZSCHILD_SC)
    {
        vr = 0.0;

        if(r > 3.0*M)
            vp = exp(-0.5/(r/M-3.0)) * sqrt(M/(r*r*r));
        else
            vp = 0.0;

        u0 = 1.0/sqrt(1-2*M/r-vr*vr/(1-2*M/r)-r*r*vp*vp);
        if(mu == 0)
            return u0;
        else if(mu == 2)
            return u0*vp;
        return 0.0;
    }
    else if(sim_Metric(theSim) == SCHWARZSCHILD_KS)
    {
        vr = -2*M/(r + 2*M);

        if(r > 3.0*M)
            vp = exp(-0.5/(r/M-3.0)) * sqrt(M/(r*r*r));
        else
            vp = 0.0;

        u0 = 1.0/sqrt(1.0-2*M/r - 4*M/r*vr - (1.0+2*M/r)*vr*vr - r*r*vp*vp);
        if(mu == 0)
            return u0;
        else if(mu == 1)
            return u0*vr;
        else if(mu == 2)
            return u0*vp;
        return 0.0;
    }
}

double metric_frame_dU_du_kep(struct Metric *g, int mu, int nu, struct Sim *theSim)
{
    double M = sim_GravM(theSim);
    double r = g->x[1];

    double u0, du0, vr, dvr, vp, dvp;

    if(mu != 1)
        return 0.0;

    if(sim_Metric(theSim) == SCHWARZSCHILD_SC)
    {
        vr = 0.0;

        if(r > 3.0*M)
        {
            vp = exp(-0.5/(r/M-3.0)) * sqrt(M/(r*r*r));
            dvp = -vp * (27*M*M-19*M*r+3*r*r) / (2*(r-3*M)*(r-3*M));
        }
        else
        {
            vp = 0.0;
            dvp = 0.0;
        }

        u0 = 1.0/sqrt(1-2*M/r-r*r*vp*vp);
        du0 = -0.5*u0*u0*u0*(2*M/(r*r) - 2*r*vp*vp - 2*r*r*vp*dvp);

        if(nu == 0)
            return du0;
        else if(nu == 2)
            return du0*vp + u0*dvp;
        return 0.0;
    }
    else if(sim_Metric(theSim) == SCHWARZSCHILD_KS)
    {
        vr = -2*M/(r + 2*M);
        dvr = 2*M/((r+2*M)*(r+2*M));

        if(r > 3.0*M)
        {
            vp = exp(-0.5/(r/M-3.0)) * sqrt(M/(r*r*r));
            dvp = -vp * (27*M*M-19*M*r+3*r*r) / (2*(r-3*M)*(r-3*M));
        }
        else
        {
            vp = 0.0;
            dvp = 0.0;
        }

        u0 = 1.0/sqrt(1.0-2*M/r - 4*M/r*vr - (1.0+2*M/r)*vr*vr - r*r*vp*vp);
        du0 = -0.5*u0*u0*u0*(2*M/(r*r) + 4*M*vr/(r*r) - 4*M/r*dvr + 2*M*vr*vr/(r*r) - 4*M/r*vr*dvr - 2*r*vp*vp - 2*r*r*vp*dvp);

        if(nu == 0)
            return du0;
        else if(nu == 1)
            return du0*vr + u0*dvr;
        else if(nu == 2)
            return du0*vp + u0*dvp;
        return 0.0;
    }
}
