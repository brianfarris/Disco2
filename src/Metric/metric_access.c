#define METRIC_PRIVATE_DEFS
#include "../Headers/Metric.h"

double metric_lapse(struct Metric *g)
{
    return g->lapse;
}

double metric_shift_u(struct Metric *g, int i)
{
    return g->shift[i];
}

double metric_shift_d(struct Metric *g, int i)
{
    return    g->gamma_dd[3*i]*g->shift[0]
            + g->gamma_uu[3*i+1]*g->shift[1] 
            + g->gamma_uu[3*i+2]*g->shift[2]; 
}

double metric_gamma_dd(struct Metric *g, int i, int j)
{
    return 0;
}

double metric_gamma_uu(struct Metric *g, int i, int j)
{
    return 0;
}

double metric_g_dd(struct Metric *g, int i, int j)
{
    return 0;
}

double metric_g_uu(struct Metric *g, int i, int j)
{
    return 0;
}

double metric_sqrtgamma(struct Metric *g)
{
    return 0;
}

double metric_sqrtg(struct Metric *g)
{
    return 0;
}
