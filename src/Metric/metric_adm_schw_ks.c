#include <math.h>
#include <stdio.h>
#include "../Headers/header.h"
#include "../Headers/Metric.h"
#include "../Headers/Sim.h"

//
//Schwarzschild metric in cylindrical Kerr-Schild coordinates
//

double metric_lapse_adm_schw_ks(double t, double r, double p, double z, 
                        struct Sim *theSim)
{
    double M = sim_GravM(theSim);
    double R = sqrt(r*r + z*z);

    return 1.0/sqrt(1.0 + 2*M/R);
}

double metric_shift_adm_schw_ks(int i, double t, double r, double p, double z, 
                        struct Sim *theSim)
{
    double M = sim_GravM(theSim);
    double R = sqrt(r*r + z*z);
    double b;

    if(i == 1)
        b = 2.0*M*r / (R*R + 2*M*R);
    else if (i == 3)
        b = 2.0*M*z / (R*R + 2*M*R);
    else
        b = 0.0;
    return b;
}

double metric_spatial_adm_schw_ks(int i, int j, double t, double r, double p, 
                            double z, struct Sim *theSim)
{
    double M = sim_GravM(theSim);
    double R = sqrt(r*r + z*z);
    double g;

    if(i == 1 && j == 1)
        g = 1.0 + 2.0*M*r*r/(R*R*R);
    else if(i == 2 && j == 2)
        g = r*r;
    else if(i == 3 && j == 3)
        g = 1.0 + 2.0*M*z*z/(R*R*R);
    else if ((i == 1 && j == 3) || (i == 3 && j == 1))
        g = 2.0*M*r*z/(R*R*R);
    else
        g = 0.0;

    return g;
}

double metric_ispatial_adm_schw_ks(int i, int j, double t, double r, double p, 
                                double z, struct Sim *theSim)
{
    double M = sim_GravM(theSim);
    double R = sqrt(r*r + z*z);
    double g;

    if(i == 1 && j == 1)
        g = (R*R*R + 2.0*M*z*z) / (R*R*R + 2.0*M*R*R);
    else if(i == 2 && j == 2)
        g = 1.0/(r*r);
    else if(i == 3 && j == 3)
        g = (R*R*R + 2.0*M*r*r) / (R*R*R + 2.0*M*R*R);
    else if ((i == 1 && j == 3) || (i == 3 && j == 1))
        g = -2.0*M*r*z / (R*R*R + 2.0*M*R*R);
    else
        g = 0.0;

    return g;
}

double metric_dlapse_adm_schw_ks(int mu, double t, double r, double p, double z, 
                            struct Sim *theSim)
{
    double M = sim_GravM(theSim);
    double R = sqrt(r*r + z*z);

    double a = 1.0/sqrt(1.0+2*M/R);
    double da;

    if(mu == 1)
        da =  a*a*a * M*r / (R*R*R);
    else if(mu == 3)
        da =  a*a*a * M*z / (R*R*R);
    else
        da = 0.0;

    return da;
}

double metric_dshift_adm_schw_ks(int mu, int i, double t, double r, double p, 
                            double z, struct Sim *theSim)
{
    double M = sim_GravM(theSim);
    double R = sqrt(r*r + z*z);

    double db;

    if(mu == 1 && i == 1)
        db =  2*M * (R*R+2*M*R-2*r*r*(1.0+M/R)) / ((R*R+2*M*R)*(R*R+2*M*R));
    else if((mu == 3 && i == 1) || (mu == 1 && i == 3))
        db =  -4*M*(R+M)*r*z / (R*(R*R+2*M*R)*(R*R+2*M*R));
    else if(mu == 3 && i == 3)
        db =  2*M * (R*R+2*M*R-2*z*z*(1.0+M/R)) / ((R*R+2*M*R)*(R*R+2*M*R));
    else
        db = 0.0;

    return db;
}

double metric_dspatial_adm_schw_ks(int mu, int i, int j, double t, double r, 
                                double p, double z, struct Sim *theSim)
{
    double M = sim_GravM(theSim);
    double R = sqrt(r*r + z*z);

    double dg;

    if(mu == 1)
    {
        if(i == 1 && j == 1)
            dg = 2*M*r * (2*R*R - 3*r*r) / (R*R*R*R*R);
        else if(i == 2 && j == 2)
            dg = 2*r;
        else if(i == 3 && j == 3)
            dg = -6*M*r*z*z / (R*R*R*R*R);
        else if((i == 1 && j == 3) || (i == 3 && j == 1))
            dg = 2*M*z * (R*R - 3*r*r) / (R*R*R*R*R);
        else
            dg = 0.0;
    }
    else if (mu == 3)
    {
        if(i == 1 && j == 1)
            dg = -6*M*r*r*z / (R*R*R*R*R);
        else if(i == 3 && j == 3)
            dg = 2*M*z * (2*R*R - 3*z*z) / (R*R*R*R*R);
        else if((i == 1 && j == 3) || (i == 3 && j == 1))
            dg = 2*M*r * (R*R - 3*z*z) / (R*R*R*R*R);
        else
            dg = 0.0;
    }
    else
        dg = 0.0;

    return dg;
}

void metric_killing_adm_schw_ks(int *k)
{
    k[0] = 1;
    k[1] = 0;
    k[2] = 1;
    k[3] = 0;
}

double metric_horizon_adm_schw_ks(struct Sim *theSim)
{
    return 2.0*sim_GravM(theSim);
}
