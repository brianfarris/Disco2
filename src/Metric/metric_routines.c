#define METRIC_PRIVATE_DEFS
#include <stdio.h>
#include "../Headers/Metric.h"

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

double metric_shear_uu(struct Metric *g, double *u, double *du, double *shear)
{
    //Fills shear tensor sigma^{ij} = shear(4*i+j) given 
    // contravariant 4-velocity u and (partial der) gradient du (d_{i}u^{j} =  du[4*i+j])

    double divu = 0.0;
    double cdu, val;
    int i,j,k;

    for(i=0; i<16; i++)
        shear[i] = 0.0;

    for(i=0; i<4; i++)
    {
        for(j=0; j<4; j++)
        {
            cdu = du[4*i+j];  //d_i u^j
            for(k=0; k<4; k++)
                cdu += metric_conn(g, j,i,k)*u[k];

            for(k=0; k<4; k++)
            {
                val = cdu*(metric_g_uu(g,i,k)+u[i]*u[k]);
                shear[4*j+k] += val;
                shear[4*k+j] += val;
            }
            if(i == j)
                divu += cdu;
        }
    }

    for(i=0; i<4; i++)
        for(j=0; j<4; j++)
            shear[4*i+j] -= 2.0/3.0*divu*(metric_g_uu(g,i,j)+u[i]*u[j]);
}

