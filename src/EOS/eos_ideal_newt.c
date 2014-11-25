#define EOS_PRIVATE_DEFS
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "../Headers/header.h"
#include "../Headers/EOS.h"
#include "../Headers/Metric.h"
#include "../Headers/Sim.h"

double eos_cs2_prim_ideal_newt(double *prim, double H, struct Sim *theSim)
{
    double GAM = sim_GAMMALAW(theSim);
    return GAM*prim[PPP]/prim[RHO];
}

double eos_eps_prim_ideal_newt(double *prim, double H, struct Sim *theSim)
{
    double GAM = sim_GAMMALAW(theSim);
    return prim[PPP]/((GAM-1.0)*prim[RHO]);
}

double eos_temp_prim_ideal_newt(double *prim, double H, struct Sim *theSim)
{
    return 0.5 * eos_mp * prim[PPP]/prim[RHO];
}
