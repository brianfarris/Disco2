#define EOS_PRIVATE_DEFS
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "../Headers/header.h"
#include "../Headers/EOS.h"
#include "../Headers/Sim.h"

double eos_cs2_gammalaw(double *prim, struct Sim *theSim)
{
    double GAM = sim_EOSPar1(theSim);
    return GAM * prim[PPP] / (1.0 + GAM*prim[PPP]/(GAM-1.0));
}

double eos_eps_gammalaw(double *prim, struct Sim *theSim)
{
    double GAM = sim_EOSPar1(theSim);
    return prim[PPP] / (GAM-1.0);
}

double eos_ppp_gammalaw(double *prim, struct Sim *theSim)
{
    return prim[RHO] * prim[PPP];
}

double eos_dpppdrho_gammalaw(double *prim, struct Sim *theSim)
{
    return prim[PPP];
}

double eos_dpppdttt_gammalaw(double *prim, struct Sim *theSim)
{
    return prim[RHO];
}

double eos_depsdrho_gammalaw(double *prim, struct Sim *theSim)
{
    return 0.0;
}

double eos_depsdttt_gammalaw(double *prim, struct Sim *theSim)
{
    double GAM = sim_EOSPar1(theSim);
    return 1.0/(GAM-1.0);
}

