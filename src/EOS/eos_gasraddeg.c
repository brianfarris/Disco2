#define EOS_PRIVATE_DEFS
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "../Headers/header.h"
#include "../Headers/EOS.h"
#include "../Headers/Sim.h"

double eos_cs2_gasraddeg(double *prim, struct Sim *theSim)
{
    double rho = prim[RHO];
    double P = eos_ppp_gasraddeg(prim, theSim);
    double e = eos_eps_gasraddeg(prim, theSim);
    double dedr = eos_depsdrho_gasraddeg(prim, theSim);
    double dedt = eos_depsdttt_gasraddeg(prim, theSim);
    double dpdr = eos_dpppdrho_gasraddeg(prim, theSim);
    double dpdt = eos_dpppdttt_gasraddeg(prim, theSim);
    double h = 1.0 + e + P/rho;

    double cs2 = (dpdr*dedt - dedr*dpdt + P*dpdt/(rho*rho)) / (dedt * h);

    return cs2;
}

double eos_eps_gasraddeg(double *prim, struct Sim *theSim)
{
    double x1 = sim_EOSPar1(theSim);
    double x2 = sim_EOSPar2(theSim);
    double x3 = sim_EOSPar3(theSim);

    double gam = sim_GAMMALAW(theSim);
    double rho = prim[RHO] * eos_rho_scale;
    double T = prim[TTT] * eos_mp*eos_c*eos_c;  //ergs

    double gas = T/((gam-1.0)*eos_mp);
    double rad = 4.0*eos_sb/(eos_c*rho)*T*T*T*T;
    double deg = eos_h*eos_c/(4*eos_mp) * pow(3*rho/(8*M_PI*eos_mp),1.0/3.0);
    
    return (x1*gas + x2*rad + x3*deg) / (eos_c*eos_c);
}

double eos_ppp_gasraddeg(double *prim, struct Sim *theSim)
{
    double x1 = sim_EOSPar1(theSim);
    double x2 = sim_EOSPar2(theSim);
    double x3 = sim_EOSPar3(theSim);

    double rho = prim[RHO] * eos_rho_scale;
    double T = prim[TTT] * eos_mp*eos_c*eos_c;  //ergs

    double gas = rho*T/eos_mp;
    double rad = 4.0*eos_sb/(3.0*eos_c)*T*T*T*T; 
    double deg = 2*M_PI*eos_h*eos_c/3.0 * pow(3*rho/(8*M_PI*eos_mp),4.0/3.0);

    return (x1*gas + x2*rad + x3*deg) / (eos_rho_scale*eos_c*eos_c);
}

double eos_dpppdrho_gasraddeg(double *prim, struct Sim *theSim)
{
    double x1 = sim_EOSPar1(theSim);
    double x2 = sim_EOSPar2(theSim);
    double x3 = sim_EOSPar3(theSim);
    double rho = prim[RHO] * eos_rho_scale;
    double T = prim[TTT] * eos_mp*eos_c*eos_c;  //ergs

    double gas = T/eos_mp;
    double rad = 0.0;
    double deg = eos_h*eos_c/(3.0*eos_mp) * pow(3*rho/(8*M_PI*eos_mp),1.0/3.0);
    
    return (x1*gas + x2*rad + x3*deg) / (eos_c*eos_c);
}

double eos_dpppdttt_gasraddeg(double *prim, struct Sim *theSim)
{
    double x1 = sim_EOSPar1(theSim);
    double x2 = sim_EOSPar2(theSim);
    double x3 = sim_EOSPar3(theSim);
    double rho = prim[RHO] * eos_rho_scale;
    double T = prim[TTT] * eos_mp*eos_c*eos_c;  //ergs

    double gas = rho;
    double rad = 16.0*eos_sb/(3.0*eos_c)*T*T*T*eos_mp;
    double deg = 0.0;
    
    return (x1*gas + x2*rad + x3*deg) / eos_rho_scale;
}

double eos_depsdrho_gasraddeg(double *prim, struct Sim *theSim)
{
    double x1 = sim_EOSPar1(theSim);
    double x2 = sim_EOSPar2(theSim);
    double x3 = sim_EOSPar3(theSim);
    double rho = prim[RHO] * eos_rho_scale;
    double T = prim[TTT] * eos_mp*eos_c*eos_c;  //ergs

    double gas = 0.0;
    double rad = -4.0*eos_sb/(eos_c*rho*rho)*T*T*T*T;
    double deg = 3*eos_h*eos_c/(32*M_PI*eos_mp*eos_mp)
                    * pow(3*rho/(8*M_PI*eos_mp),-2.0/3.0);
    
    return (x1*gas + x2*rad + x3*deg) * (eos_rho_scale / (eos_c*eos_c));
}

double eos_depsdttt_gasraddeg(double *prim, struct Sim *theSim)
{
    double x1 = sim_EOSPar1(theSim);
    double x2 = sim_EOSPar2(theSim);
    double x3 = sim_EOSPar3(theSim);

    double gam = sim_GAMMALAW(theSim);
    double rho = prim[RHO] * eos_rho_scale;
    double T = prim[TTT] * eos_mp*eos_c*eos_c;  //ergs

    double gas = 1.0/(gam-1.0);
    double rad = 16.0*eos_sb/(eos_c*rho)*T*T*T*eos_mp;
    double deg = 0.0;
    
    return (x1*gas + x2*rad + x3*deg);
}

