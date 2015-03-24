#include <stdio.h>
#include "../Headers/header.h"
#include "../Headers/Metric.h"
#include "../Headers/Sim.h"

//
//Flat space-time metric in cylindrical coordinates
//

double metric_lapse_adm_sr(double t, double r, double p, double z, 
                        struct Sim *theSim)
{
    return 1.0;
}

double metric_shift_adm_sr(int i, double t, double r, double p, double z, 
                        struct Sim *theSim)
{
    return 0.0;
}

double metric_spatial_adm_sr(int i, int j, double t, double r, double p, 
                            double z, struct Sim *theSim)
{
    double g = 0.0;
    if(i != j)
        g = 0.0;
    else if (i == 1 || i == 3)
        g = 1.0;
    else if (i == 2)
        g = r*r;
    return g;
}

double metric_ispatial_adm_sr(int i, int j, double t, double r, double p, 
                                double z, struct Sim *theSim)
{
    double g = 0.0;
    if(i != j)
        g = 0.0;
    else if (i == 1 || i == 3)
        g = 1.0;
    else if (i == 2)
        g = 1.0/(r*r);
    return g;
}

double metric_dlapse_adm_sr(int mu, double t, double r, double p, double z, 
                            struct Sim *theSim)
{
    if(PRINTTOOMUCH)
        printf("in ADM SR dlapse\n");
    return 0.0;
}

double metric_dshift_adm_sr(int mu, int i, double t, double r, double p, 
                            double z, struct Sim *theSim)
{
    return 0.0;
}

double metric_dspatial_adm_sr(int mu, int i, int j, double t, double r, 
                                double p, double z, struct Sim *theSim)
{
    double dg;

    if(mu == 1 && i == 2 && j == 2)
        dg = 2*r;
    else
        dg = 0.0;

    return dg;
}

void metric_killing_adm_sr(int *k)
{
    k[0] = 1;
    k[1] = 0;
    k[2] = 1;
    k[3] = 1;
}

double metric_horizon_adm_sr(struct Sim *theSim)
{
    return -1.0;
}
