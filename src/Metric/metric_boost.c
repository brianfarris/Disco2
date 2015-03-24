#define METRIC_PRIVATE_DEFS
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "../Headers/header.h"
#include "../Headers/Sim.h"
#include "../Headers/Metric.h"

double metric_boost_none(int i, double t, double r, double p, double z, 
                        struct Sim *theSim)
{
    return 0;
}

double metric_dboost_none(int mu, int i, double t, double r, double p, 
                        double z, struct Sim *theSim)
{
    return 0;
}

void metric_killing_boost_none(int *k){}

double metric_boost_rigid(int i, double t, double r, double p, double z, 
                            struct Sim *theSim)
{
    double b;
    if(i == 2)
        b = sim_BinW(theSim);
    else
        b = 0.0;
    return b;
}

double metric_dboost_rigid(int mu, int i, double t, double r, double p, 
                        double z, struct Sim *theSim)
{
    return 0;
}

void metric_killing_boost_rigid(int *k){}

double metric_boost_bin(int i, double t, double r, double p, double z, 
                            struct Sim *theSim)
{
    double a = sim_BinA(theSim);
    double w = sim_BinW(theSim);

    double b;
    if(i == 1)
        b = a*w*sin(p);
    else if (i == 2)
        b = w + a*w*cos(p)/r;
    else
        b = 0.0;

    return b;
}

double metric_dboost_bin(int mu, int i, double t, double r, double p, 
                        double z, struct Sim *theSim)
{
    double a = sim_BinA(theSim);
    double w = sim_BinW(theSim);

    double db;

    if(mu == 1 && i == 2)
        db = -a*w*cos(p)/(r*r);
    else if (mu == 2 && i == 1)
        db = a*w*cos(p);
    else if (mu == 2 && i == 2)
        db = -a*w*sin(p)/r;
    else
        db = 0;

    return db;
}

void metric_killing_boost_bin(int *k)
{
    k[2] = 0;
}
