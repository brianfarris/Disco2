#include "../Headers/Metric.h"

//Flat space-time metric in cartesian coordinates.
//
double metric_g_dd_exact_sr_cart(int mu, int nu, double t, double r, double p, double z, struct Sim *theSim)
{
    if(mu != nu)
        return 0.0;
    if(mu == 0)
        return -1.0;
    if(mu == 1 || mu == 2 || mu == 3)
        return 1.0;
    return 0.0;
}

double metric_g_uu_exact_sr_cart(int mu, int nu, double t, double r, double p, double z, struct Sim *theSim)
{
    if(mu != nu)
        return 0.0;
    if(mu == 0)
        return -1.0;
    if(mu == 1 || mu == 2 || mu == 3)
        return 1.0;
    return 0.0;
}

double metric_dg_dd_exact_sr_cart(int k, int mu, int nu, double t, double r, double p, double z, struct Sim *theSim)
{
    return 0.0;
}

double metric_dg_uu_exact_sr_cart(int k, int mu, int nu, double t, double r, double p, double z, struct Sim *theSim)
{
    return 0.0;
}

void metric_killing_exact_sr_cart(int *k)
{
    k[0] = 1;
    k[1] = 1;
    k[2] = 1;
    k[3] = 1;
}

double metric_horizon_exact_sr_cart(struct Sim *theSim)
{
    return -1.0;
}
