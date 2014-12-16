#define EOS_PRIVATE_DEFS
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "../Headers/header.h"
#include "../Headers/EOS.h"
#include "../Headers/Sim.h"

double eos_cs2_gasrad(double *prim, struct Sim *theSim)
{
    double rho = prim[RHO];
    double T = prim[TTT];
    double P = eos_ppp_gasrad(prim, theSim);
    double e = eos_eps_gasrad(prim, theSim);
    double dedt = eos_depsdttt_gasrad(prim, theSim);
    double dpdr = eos_dpppdrho_gasrad(prim, theSim);
    double dpdt = eos_dpppdttt_gasrad(prim, theSim);
    double h = 1.0 + e + P/rho;

    double cs2 = (dpdr*dedt + T*dpdt*dpdt/(rho*rho)) / (dedt * h);

    return cs2;
}

double eos_eps_gasrad(double *prim, struct Sim *theSim)
{
    double x1 = sim_EOSPar1(theSim);
    double x2 = sim_EOSPar1(theSim);
    double gam = sim_GAMMALAW(theSim);
    double rho = prim[RHO] * eos_rho_scale;
    double T = prim[TTT] * eos_mp*eos_c*eos_c;  //ergs

    double gas = T/((gam-1.0)*eos_mp);
    double rad = 4.0*eos_sb/(eos_c*rho)*T*T*T*T;
    
    return (x1*gas + x2*rad) / (eos_c*eos_c);
}

double eos_ppp_gasrad(double *prim, struct Sim *theSim)
{
    double x1 = sim_EOSPar1(theSim);
    double x2 = sim_EOSPar1(theSim);
    double rho = prim[RHO] * eos_rho_scale;
    double T = prim[TTT] * eos_mp*eos_c*eos_c;  //ergs

    double gas = rho*T/eos_mp;
    double rad = 4.0*eos_sb/(3.0*eos_c)*T*T*T*T;
    
    return (x1*gas + x2*rad) / (eos_rho_scale*eos_c*eos_c);
}

double eos_dpppdrho_gasrad(double *prim, struct Sim *theSim)
{
    double x1 = sim_EOSPar1(theSim);
    double x2 = sim_EOSPar1(theSim);
    double T = prim[TTT] * eos_mp*eos_c*eos_c;  //ergs

    double gas = T/eos_mp;
    double rad = 0.0;
    
    return (x1*gas + x2*rad) / (eos_c*eos_c);
}

double eos_dpppdttt_gasrad(double *prim, struct Sim *theSim)
{
    double x1 = sim_EOSPar1(theSim);
    double x2 = sim_EOSPar1(theSim);
    double rho = prim[RHO] * eos_rho_scale;
    double T = prim[TTT] * eos_mp*eos_c*eos_c;  //ergs

    double gas = rho;
    double rad = 16.0*eos_sb/(3.0*eos_c)*T*T*T*eos_mp;
    
    return (x1*gas + x2*rad) / eos_rho_scale;
}

double eos_depsdrho_gasrad(double *prim, struct Sim *theSim)
{
    double x1 = sim_EOSPar1(theSim);
    double x2 = sim_EOSPar1(theSim);
    double rho = prim[RHO] * eos_rho_scale;
    double T = prim[TTT] * eos_mp*eos_c*eos_c;  //ergs

    double gas = 0.0;
    double rad = -4.0*eos_sb/(eos_c*rho*rho)*T*T*T*T;
    
    return (x1*gas + x2*rad) * (eos_rho_scale / (eos_c*eos_c));
}

double eos_depsdttt_gasrad(double *prim, struct Sim *theSim)
{
    double x1 = sim_EOSPar1(theSim);
    double x2 = sim_EOSPar1(theSim);
    double gam = sim_GAMMALAW(theSim);
    double rho = prim[RHO] * eos_rho_scale;
    double T = prim[TTT] * eos_mp*eos_c*eos_c;  //ergs

    double gas = 1.0/(gam-1.0);
    double rad = 16.0*eos_sb/(eos_c*rho)*T*T*T*eos_mp;
    
    return (x1*gas + x2*rad);
}

