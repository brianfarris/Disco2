#include "../Headers/Metric.h"

//Flat space-time metric in cylindrical coordinates
//
double metric_g_dd_exact_sr(int mu, int nu, double t, double r, double p, double z, struct Sim *theSim)
{
    if(mu != nu)
        return 0.0;
    if(mu == 0)
        return -1.0;
    if(mu == 1)
        return 1.0;
    if(mu == 2)
        return r*r;
    if(mu == 3)
        return 1.0;
    return 0.0;
}

double metric_g_uu_exact_sr(int mu, int nu, double t, double r, double p, double z, struct Sim *theSim)
{
    if(mu != nu)
        return 0.0;
    if(mu == 0)
        return -1.0;
    if(mu == 1)
        return 1.0;
    if(mu == 2)
        return 1.0/(r*r);
    if(mu == 3)
        return 1.0;
    return 0.0;
}

double metric_dg_dd_exact_sr(int k, int mu, int nu, double t, double r, double p, double z, struct Sim *theSim)
{
    if(k==1 && mu==2 && nu==2)
        return 2.0*r;
    return 0.0;
}

double metric_dg_uu_exact_sr(int k, int mu, int nu, double t, double r, double p, double z, struct Sim *theSim)
{
    if(k==1 && mu==2 && nu==2)
        return -2.0/(r*r*r);
    return 0.0;
}

void metric_killing_exact_sr(int *k)
{
    k[0] = 1;
    k[1] = 0;
    k[2] = 1;
    k[3] = 1;
}

double metric_horizon_exact_sr(struct Sim *theSim)
{
    return -1.0;
}
