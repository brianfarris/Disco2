#define METRIC_PRIVATE_DEFS
#include <math.h>
#include "../Headers/Metric.h"
#include "../Headers/header.h"
#include "../Headers/Sim.h"

double metric_lapse(struct Metric *g)
{
    return sqrt(-1.0/(g->g_uu[0]));
}

double metric_shift_u(struct Metric *g, int i)
{
    return -(g->g_uu[i+1])/(g->g_uu[0]);
}

double metric_shift_d(struct Metric *g, int i)
{
    return g->g_dd[i+1];
}

double metric_gamma_dd(struct Metric *g, int i, int j)
{
    return metric_g_dd(g,i+1,j+1);
}

double metric_gamma_uu(struct Metric *g, int i, int j)
{
    return metric_g_uu(g,i+1,j+1) - g->g_uu[i+1]*g->g_uu[j+1]/(g->g_uu[0]);
}

double metric_g_dd(struct Metric *g, int i, int j)
{
    if(i <= j)
        return g->g_dd[j+4*i-i*(i+1)/2];
    return g->g_dd[i+4*j-j*(j+1)/2];
}

double metric_g_uu(struct Metric *g, int i, int j)
{
    if(i <= j)
        return g->g_uu[j+4*i-i*(i+1)/2];
    return g->g_uu[i+4*j-j*(j+1)/2];
}

double metric_sqrtgamma(struct Metric *g)
{
    return sqrt(g->g_dd[4]*(g->g_dd[7]*g->g_dd[9] - g->g_dd[8]*g->g_dd[8])
        + g->g_dd[5]*(g->g_dd[8]*g->g_dd[6] - g->g_dd[5]*g->g_dd[9])
        + g->g_dd[6]*(g->g_dd[5]*g->g_dd[8] - g->g_dd[7]*g->g_dd[6]));
}

double metric_sqrtg(struct Metric *g)
{
    double det = 0.0;

    det = g->g_dd[0] * metric_sqrtgamma(g);
    if(g->g_dd[1] != 0.0)
        det -= g->g_dd[1]*(g->g_dd[1]*(g->g_dd[7]*g->g_dd[9] - g->g_dd[8]*g->g_dd[8])
                        + g->g_dd[5]*(g->g_dd[8]*g->g_dd[3] - g->g_dd[2]*g->g_dd[9])
                        + g->g_dd[6]*(g->g_dd[2]*g->g_dd[8] - g->g_dd[7]*g->g_dd[3]));
    if(g->g_dd[2] != 0.0)
        det += g->g_dd[2]*(g->g_dd[1]*(g->g_dd[5]*g->g_dd[9] - g->g_dd[8]*g->g_dd[6])
                        + g->g_dd[4]*(g->g_dd[8]*g->g_dd[3] - g->g_dd[2]*g->g_dd[9])
                        + g->g_dd[6]*(g->g_dd[2]*g->g_dd[6] - g->g_dd[5]*g->g_dd[3]));
    if(g->g_dd[3] != 0.0)
        det -= g->g_dd[3]*(g->g_dd[1]*(g->g_dd[5]*g->g_dd[8] - g->g_dd[7]*g->g_dd[6])
                        + g->g_dd[4]*(g->g_dd[7]*g->g_dd[3] - g->g_dd[2]*g->g_dd[8])
                        + g->g_dd[5]*(g->g_dd[2]*g->g_dd[6] - g->g_dd[5]*g->g_dd[3]));
    return sqrt(-det);
}

double metric_dg_dd(struct Metric *g, int k, int i, int j)
{
    if(g->killing[k])
        return 0.0;
    if(g->length_dg == 0)
        metric_create_der(g);
    int a, b;
    a = 0;
    for(b = 0; b<k; b++)
        if(!(g->killing[b]))
            a++;
    if(i <= j)
        return g->dg_dd[10*a + j+4*i-i*(i+1)/2];
    return g->dg_dd[10*a + i+4*j-j*(j+1)/2];
}

double metric_dg_uu(struct Metric *g, int k, int i, int j)
{
    if(g->killing[k])
        return 0.0;
    if(g->length_dg == 0)
        metric_create_der(g);
    int a, b;
    a = 0;
    for(b = 0; b<k; b++)
        if(!g->killing[b])
            a++;
    if(i <= j)
        return g->dg_uu[10*a + j+4*i-i*(i+1)/2];
    return g->dg_uu[10*a + i+4*j-j*(j+1)/2];
}

double metric_dlapse(struct Metric *g, int k)
{
    return 0.5 * metric_dg_uu(g,k,0,0) / sqrt(-g->g_uu[0]*g->g_uu[0]*g->g_uu[0]);
}

double metric_dshift_u(struct Metric *g, int k, int i)
{
    if(i < 0 || i > 2)
        return 0.0;
    return g->g_uu[i+1]*metric_dg_uu(g,k,0,0)/(g->g_uu[0]*g->g_uu[0]) - metric_dg_uu(g,k,0,i+1)/(g->g_uu[0]);
}

int metric_killcoord(struct Metric *g, int k)
{
    return g->killing[k];
}

double metric_conn(struct Metric *g, int tau, int mu, int nu)
{
    //Returns value of the connection coefficient \Gamma^\tau_{\mu\nu}
    //
    //If you need to calculate LOTS of christoffel symbols, this is not the
    //most efficient way.
    
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

    //Calculate 4-velocity u.
    a = metric_lapse(g);
    for(mu=0; mu<3; mu++)
        b[mu] = metric_shift_u(g, mu);

    u[0] = 1.0/sqrt(-metric_g_dd(g,0,0) - 2*metric_dot3_u(g,b,v)-metric_square3_u(g,v));
    u[1] = u[0]*v[0];
    u[2] = u[0]*v[1];
    u[3] = u[0]*v[2];

    //Calculate coordinate derivatives of 4-velocity.
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

    
    //Calculate all Christoffel Symbols
    double christoffel_d[64];
    for(mu=0; mu<64; mu++)
        christoffel_d[mu] = 0.0;

    double christoffel[64];
    for(mu=0; mu<64; mu++)
        christoffel[mu] = 0.0;

    //First calculate fully lowered symbols.
    for(la=0; la<4; la++)
        for(mu=0; mu<4; mu++)
        {
            for(nu=0; nu<mu; nu++)
            {
                double dg = metric_dg_dd(g, la, mu,nu);
                if(dg != 0.0)
                {
                    christoffel_d[16*la+4*mu+nu] -= dg;
                    christoffel_d[16*la+4*nu+mu] -= dg;
                    christoffel_d[16*mu+4*la+nu] += dg;
                    christoffel_d[16*mu+4*nu+la] += dg;
                    christoffel_d[16*nu+4*mu+la] += dg;
                    christoffel_d[16*nu+4*la+mu] += dg;
                }
            }
            double dg = metric_dg_dd(g, la, mu,mu);
            if(dg != 0.0)
            {
                christoffel_d[16*la+4*mu+mu] -= dg;
                christoffel_d[16*mu+4*la+mu] += dg;
                christoffel_d[16*mu+4*mu+la] += dg;
            }
        }
   
    //Now raise the first index to get the standard symbols.
    for(mu=0; mu<4; mu++)
    {
        for(nu=0; nu<mu; nu++)
        {
            double ig = metric_g_uu(g, mu, nu);

            if(ig != 0.0)
            {
                int si;
                for(la=0; la<4; la++)
                    for(si=0; si<4; si++)
                        christoffel[16*mu+4*la+si] += ig*christoffel_d[16*nu+4*la+si];
                for(la=0; la<4; la++)
                    for(si=0; si<4; si++)
                        christoffel[16*nu+4*la+si] += ig*christoffel_d[16*mu+4*la+si];
            }
        }
        double ig = metric_g_uu(g, mu, mu);

        if(ig != 0.0)
        {
            int si;
            for(la=0; la<4; la++)
                for(si=0; si<4; si++)
                    christoffel[16*mu+4*la+si] += ig*christoffel_d[16*mu+4*la+si];
        }
    }

    for(mu=0; mu<64; mu++)
        christoffel[mu] *= 0.5;

    //Add transverse components to shear
    for(mu=0; mu<4; mu++)
    {
        for(nu=0; nu<4; nu++)
        {
            cdu = du[4*mu+nu];  //d_mu u^nu
            for(la=0; la<4; la++)
            {
             //   cdu += metric_conn(g, nu,mu,la)*u[la];
                cdu += christoffel[16*nu+4*mu+la]*u[la];
            }

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

    //Subtract trace to make traceless!
    if(sim_N_global(theSim, Z_DIR) != 1)
        divu *= 2.0/3.0;

    for(mu=0; mu<4; mu++)
        for(nu=0; nu<4; nu++)
            shear[4*mu+nu] -= divu*(metric_g_uu(g,mu,nu)+u[mu]*u[nu]);
}
