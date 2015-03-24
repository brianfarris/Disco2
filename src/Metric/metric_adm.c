#define METRIC_PRIVATE_DEFS
#include <stdio.h>
#include <stdlib.h>
#include "../Headers/header.h"
#include "../Headers/Sim.h"
#include "../Headers/Metric.h"

double metric_shift_adm(int i, double t, double r, double p, double z, 
                        struct Sim *theSim)
{
    return metric_shift_adm_exact(i,t,r,p,z,theSim) + 
            metric_shift_adm_boost(i,t,r,p,z,theSim);
}

double metric_dshift_adm(int mu, int i, double t, double r, double p, 
                        double z, struct Sim *theSim)
{
    return metric_dshift_adm_exact(mu,i,t,r,p,z,theSim) + 
            metric_dshift_adm_boost(mu,i,t,r,p,z,theSim);
}
