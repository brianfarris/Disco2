#define EOS_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include "../Headers/header.h"
#include "../Headers/EOS.h"
#include "../Headers/Metric.h"

double eos_cool_none(double *x, double *prim, struct Sim *theSim)
{
    return 0.0;
}

double eos_cool_isotherm(double *x, double *prim, struct Sim *theSim)
{
    double T0 = sim_CoolPar1(theSim);
    double tau = sim_CoolPar2(theSim);

    return (prim[PPP]/prim[RHO] - T0) / tau;
}

double eos_cool_bb_es(double *x, double *prim, struct Sim *theSim)
{
    //TODO: Work in 2D and 3D?
    double temp = prim[PPP] / prim[RHO];
    double q0 = sim_CoolPar1(theSim);
    return q0 * temp*temp*temp*temp / prim[RHO];
}

double eos_cool_visc(double *x, double *prim, struct Sim *theSim)
{
    return 0.0;
}
