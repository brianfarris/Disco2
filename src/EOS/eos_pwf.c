#define EOS_PRIVATE_DEFS
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "../Headers/header.h"
#include "../Headers/EOS.h"
#include "../Headers/Sim.h"

// Equation of state which matches Popham, Woosley, and Fryer 1999 (PWF).
// The gas is composed of free nucleons, helium, relativistic elections & 
// positrons, and photons.  The nucleons and helium contribute gas pressure,
// the photons, electrons, and positrons contribute radiation pressure, and 
// the electrons contribute relativistic degeneracy pressure.
//
// The explicit expressions for reletivistic degneracy pressure and radiation
// are used.  This is only a good approximation if the temperature is greater
// than twice the electron mass, so e+e- pairs can be spontaneously generated.
// Chemical and nuclear equilibria are not calculated explicity, rather the
// mass fraction of nucleons Xnuc is calculated according to equation 3.10 in
// PWF.
//
// This setup is not meant to be an accurate depiction of the real world, but 
// rather as a vastly simplified model which can be directly compared to PWF
// for benchmarking.

double eos_eps_pwf_loc(double rho, double T, double xnuc, double gam, 
                        double x1, double x2, double x3);
double eos_ppp_pwf_loc(double rho, double T, double xnuc, double gam, 
                        double x1, double x2, double x3);
double eos_dpppdrho_pwf_loc(double rho, double T, double xnuc, double gam, 
                        double x1, double x2, double x3);
double eos_dpppdttt_pwf_loc(double rho, double T, double xnuc, double gam, 
                        double x1, double x2, double x3);
double eos_depsdrho_pwf_loc(double rho, double T, double xnuc, double gam, 
                        double x1, double x2, double x3);
double eos_depsdttt_pwf_loc(double rho, double T, double xnuc, double gam, 
                        double x1, double x2, double x3);

double eos_cs2_pwf(double *prim, struct Sim *theSim)
{
    double x1 = sim_EOSPar1(theSim);
    double x2 = sim_EOSPar2(theSim);
    double x3 = sim_EOSPar3(theSim);

    double gam = sim_GAMMALAW(theSim);
    double rho = prim[RHO];
    double rhocgs = prim[RHO] * eos_rho_scale;
    double Tcgs = prim[TTT] * eos_mp*eos_c*eos_c;  //ergs
    double xnuc = eos_xnuc_pwf(rhocgs, Tcgs/eos_k);

    double P = eos_ppp_pwf_loc(rhocgs, Tcgs, xnuc, gam, x1, x2, x3);
    double e = eos_eps_pwf_loc(rhocgs, Tcgs, xnuc, gam, x1, x2, x3);
    double dedr = eos_depsdrho_pwf_loc(rhocgs, Tcgs, xnuc, gam, x1, x2, x3);
    double dedt = eos_depsdttt_pwf_loc(rhocgs, Tcgs, xnuc, gam, x1, x2, x3);
    double dpdr = eos_dpppdrho_pwf_loc(rhocgs, Tcgs, xnuc, gam, x1, x2, x3);
    double dpdt = eos_dpppdttt_pwf_loc(rhocgs, Tcgs, xnuc, gam, x1, x2, x3);
    double h = 1.0 + e + P/rho;

    double cs2 = (dpdr*dedt - dedr*dpdt + P*dpdt/(rho*rho)) / (dedt * h);

    return cs2;
}

double eos_eps_pwf(double *prim, struct Sim *theSim)
{
    double x1 = sim_EOSPar1(theSim);
    double x2 = sim_EOSPar2(theSim);
    double x3 = sim_EOSPar3(theSim);

    double gam = sim_GAMMALAW(theSim);
    double rho = prim[RHO] * eos_rho_scale;
    double T = prim[TTT] * eos_mp*eos_c*eos_c;  //ergs
    double xnuc = eos_xnuc_pwf(rho, T/eos_k);

    return eos_eps_pwf_loc(rho, T, xnuc, gam, x1, x2, x3);
}

double eos_ppp_pwf(double *prim, struct Sim *theSim)
{
    double x1 = sim_EOSPar1(theSim);
    double x2 = sim_EOSPar2(theSim);
    double x3 = sim_EOSPar3(theSim);

    double gam = sim_GAMMALAW(theSim);
    double rho = prim[RHO] * eos_rho_scale;
    double T = prim[TTT] * eos_mp*eos_c*eos_c;  //ergs
    double xnuc = eos_xnuc_pwf(rho, T/eos_k);

    return eos_ppp_pwf_loc(rho, T, xnuc, gam, x1, x2, x3);
}

double eos_dpppdrho_pwf(double *prim, struct Sim *theSim)
{
    double x1 = sim_EOSPar1(theSim);
    double x2 = sim_EOSPar2(theSim);
    double x3 = sim_EOSPar3(theSim);

    double gam = sim_GAMMALAW(theSim);
    double rho = prim[RHO] * eos_rho_scale;
    double T = prim[TTT] * eos_mp*eos_c*eos_c;  //ergs
    double xnuc = eos_xnuc_pwf(rho, T/eos_k);

    return eos_dpppdrho_pwf_loc(rho, T, xnuc, gam, x1, x2, x3);
}

double eos_dpppdttt_pwf(double *prim, struct Sim *theSim)
{
    double x1 = sim_EOSPar1(theSim);
    double x2 = sim_EOSPar2(theSim);
    double x3 = sim_EOSPar3(theSim);

    double gam = sim_GAMMALAW(theSim);
    double rho = prim[RHO] * eos_rho_scale;
    double T = prim[TTT] * eos_mp*eos_c*eos_c;  // (ergs)
    double xnuc = eos_xnuc_pwf(rho, T/eos_k);

    return eos_dpppdttt_pwf_loc(rho, T, xnuc, gam, x1, x2, x3);
}

double eos_depsdrho_pwf(double *prim, struct Sim *theSim)
{
    double x1 = sim_EOSPar1(theSim);
    double x2 = sim_EOSPar2(theSim);
    double x3 = sim_EOSPar3(theSim);
    
    double gam = sim_GAMMALAW(theSim);
    double rho = prim[RHO] * eos_rho_scale;
    double T = prim[TTT] * eos_mp*eos_c*eos_c;  //ergs
    double xnuc = eos_xnuc_pwf(rho, T/eos_k);

    return eos_depsdrho_pwf_loc(rho, T, xnuc, gam, x1, x2, x3);
}

double eos_depsdttt_pwf(double *prim, struct Sim *theSim)
{
    double x1 = sim_EOSPar1(theSim);
    double x2 = sim_EOSPar2(theSim);
    double x3 = sim_EOSPar3(theSim);

    double gam = sim_GAMMALAW(theSim);
    double rho = prim[RHO] * eos_rho_scale;
    double T = prim[TTT] * eos_mp*eos_c*eos_c;  //ergs
    double xnuc = eos_xnuc_pwf(rho, T/eos_k);

    return eos_depsdttt_pwf_loc(rho, T, xnuc, gam, x1, x2, x3);
}

double eos_xnuc_pwf(double rhocgs, double Tk)
{
// Returns the mass fraction of nucleons X_nuc, calculated from eq 3.10 in 
// PWF.  Requires rho in g/cm^3 and T in Kelvin.
    double xnuc = 30.97 * pow(rhocgs*1.0e-10,-0.75) * pow(Tk*1.0e-11,1.125)
                    * exp(-6.096/(Tk*1.0e-11));
    if(xnuc > 1.0)
        return 1.0;
    return xnuc;
}

double eos_dxnudrho_pwf(double rhocgs, double Tk, double xnuc)
{
    //Returns d(X_nuc)/dRho in CGS units.

    if(xnuc > 1.0)
        return 0.0;

    return -0.75 * xnuc / rhocgs;
}

double eos_dxnudttt_pwf(double rhocgs, double Tk, double xnuc)
{
    //Returns d(X_nuc)/dT in CGS units.
    if(xnuc > 1.0)
        return 0.0;

    return (1.125/Tk + 6.096/(Tk*Tk*1.0e-11))*xnuc;
}

//
// Local functions, not meant to be used outside this file.
//

double eos_eps_pwf_loc(double rho, double T, double xnuc, double gam, 
                        double x1, double x2, double x3)
{
    double gas = (0.25+0.75*xnuc)/((gam-1.0)*eos_mp) * T;
    double rad = (11.0/4.0) * 4.0*eos_sb/(eos_c*rho) * T*T*T*T;
    double deg = 3*eos_h*eos_c/(4*eos_mp)*pow(3*rho/(16*M_PI*eos_mp),1.0/3.0);
    
    return (x1*gas + x2*rad + x3*deg) / (eos_c*eos_c);
}

double eos_ppp_pwf_loc(double rho, double T, double xnuc, double gam, 
                        double x1, double x2, double x3)
{
    double gas = (0.25+0.75*xnuc)*rho*T/eos_mp;
    double rad = (11.0/12.0) * 4.0*eos_sb/(eos_c) * T*T*T*T; 
    double deg = 2*M_PI*eos_h*eos_c/3.0 * pow(3*rho/(16*M_PI*eos_mp),4.0/3.0);

    return (x1*gas + x2*rad + x3*deg) / (eos_rho_scale*eos_c*eos_c);
}

double eos_dpppdrho_pwf_loc(double rho, double T, double xnuc, double gam, 
                        double x1, double x2, double x3)
{
    double dxdr = eos_dxnudrho_pwf(rho, T/eos_k, xnuc);

    double gas = 0.75*dxdr*rho*T/eos_mp + (0.25+0.75*xnuc)*T/eos_mp;
    double rad = 0.0;
    double deg = eos_h*eos_c/(6.0*eos_mp)
                    * pow(3*rho/(16*M_PI*eos_mp),1.0/3.0);
    
    return (x1*gas + x2*rad + x3*deg) / (eos_c*eos_c);
}

double eos_dpppdttt_pwf_loc(double rho, double T, double xnuc, double gam, 
                        double x1, double x2, double x3)
{
    double dxdt = eos_dxnudttt_pwf(rho, T/eos_k, xnuc)/eos_k; // (#/ergs)

    double gas = 0.75*dxdt*rho*T/eos_mp + (0.25+0.75*xnuc)*rho/eos_mp;
    double rad = (11.0/3.0) * 4.0*eos_sb/eos_c * T*T*T;
    double deg = 0.0;
    
    return (x1*gas + x2*rad + x3*deg) * (eos_mp/eos_rho_scale);
}

double eos_depsdrho_pwf_loc(double rho, double T, double xnuc, double gam, 
                        double x1, double x2, double x3)
{
    double dxdr = eos_dxnudrho_pwf(rho, T/eos_k, xnuc); // (cm^3/g)

    double gas = 0.75*dxdr / ((gam-1.0)*eos_mp) * T;
    double rad = (-1.0) * 4.0*eos_sb/(eos_c*rho*rho) * T*T*T*T;
    double deg = 3.0*eos_h*eos_c/(64*M_PI*eos_mp*eos_mp)
                    * pow(3*rho/(16*M_PI*eos_mp),-2.0/3.0);
    
    return (x1*gas + x2*rad + x3*deg) * (eos_rho_scale / (eos_c*eos_c));
}

double eos_depsdttt_pwf_loc(double rho, double T, double xnuc, double gam, 
                        double x1, double x2, double x3)
{
    double dxdt = eos_dxnudttt_pwf(rho, T/eos_k, xnuc)/eos_k; // (#/erg)

    double gas = (0.75*dxdt*T + (0.25+0.75*xnuc))/((gam-1.0)*eos_mp);
    double rad = (11.0) * 4.0*eos_sb/(eos_c*rho)*T*T*T;
    double deg = 0.0;
    
    return (x1*gas + x2*rad + x3*deg) * eos_mp;
}
