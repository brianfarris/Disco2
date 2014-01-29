#define METRIC_PRIVATE_DEFS
#include <math.h>
#include "../Headers/Metric.h"

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
    if(i <= 0 || i > 3)
        return 0.0;
    return g->g_uu[i]*metric_dg_uu(g,k,0,0)/(g->g_uu[0]*g->g_uu[0]) - metric_dg_uu(g,k,0,i)/(g->g_uu[0]);
}

int metric_killcoord(struct Metric *g, int k)
{
    return g->killing[k];
}
