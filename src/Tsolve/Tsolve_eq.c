#define TSOLVE_PRIVATE_DEFS
#include "../Headers/Tsolve.h"
#include <math.h>
double tsolve_E_eq (double T, void *params)
{
    struct Tsolve *theTsolve = (struct Tsolve *) params;

    double rhoe = theTsolve->rhoe_or_Ptot;
    double Sigma = theTsolve->Sigma;
    double Gamma = theTsolve->Gamma;
    double xi_r = theTsolve->xi_r;
    double Omega_eff = theTsolve->Omega_eff;
    double n = 1./(Gamma-1.);
    double absT = fabs(T);
    double cs = sqrt(1./3.*rhoe/Sigma+(1.-n/3.)*absT);
    double h = cs/Omega_eff;
    return rhoe - n*Sigma*absT - 3.*xi_r* h *pow(absT,4);
}

double tsolve_E_eq_deriv (double T, void *params)
{
    struct Tsolve *theTsolve = (struct Tsolve *) params;

    double rhoe = theTsolve->rhoe_or_Ptot;
    double Sigma = theTsolve->Sigma;
    double Gamma = theTsolve->Gamma;
    double xi_r = theTsolve->xi_r;
    double Omega_eff = theTsolve->Omega_eff;

    double n = 1./(Gamma-1.);
    double cs = sqrt(1./3.*rhoe/Sigma + (1.-n/3.)*T);
    double h = cs/Omega_eff;
 
    double dcs_dT_o_cs = 0.5/cs/cs*(1.-n/3.);

    return - n*Sigma - 3.*xi_r*h*pow(T,4) * (dcs_dT_o_cs + 4./T);
}

void tsolve_E_eq_fdf (double T, void *params, 
        double *y, double *dy)
{
    struct Tsolve *theTsolve = (struct Tsolve *) params;

    double rhoe = theTsolve->rhoe_or_Ptot;
    double Sigma = theTsolve->Sigma;
    double Gamma = theTsolve->Gamma;
    double xi_r = theTsolve->xi_r;
    double Omega_eff = theTsolve->Omega_eff;

    double n = 1./(Gamma-1.);
    double cs = sqrt(1./3.*rhoe/Sigma + (1.-n/3.)*T);
    double h = cs/Omega_eff;
    double dcs_dT_o_cs = 0.5/cs/cs*(1.-n/3.);


    *y = rhoe - n*Sigma*T - 3.*xi_r* h *pow(T,4);

    *dy = - n*Sigma - 3.*xi_r*h*pow(T,4) * (dcs_dT_o_cs + 4./T);
}

double tsolve_P_eq (double T, void *params)
{
    struct Tsolve *theTsolve = (struct Tsolve *) params;

    double Ptot = theTsolve->rhoe_or_Ptot;
    double Sigma = theTsolve->Sigma;
    double xi_r = theTsolve->xi_r;
    double Omega_eff = theTsolve->Omega_eff;
    double cs = sqrt(Ptot/Sigma);
    double h = cs/Omega_eff;
    return Ptot - Sigma*T - xi_r* h *pow(T,4);
}

double tsolve_P_eq_deriv (double T, void *params)
{
    struct Tsolve *theTsolve = (struct Tsolve *) params;

    double Ptot = theTsolve->rhoe_or_Ptot;
    double Sigma = theTsolve->Sigma;
    double xi_r = theTsolve->xi_r;
    double Omega_eff = theTsolve->Omega_eff;

    double cs = sqrt(Ptot/Sigma);
    double h = cs/Omega_eff;
 

    return - Sigma - 4.*xi_r*h*pow(T,3);
}

void tsolve_P_eq_fdf (double T, void *params, 
        double *y, double *dy)
{
    struct Tsolve *theTsolve = (struct Tsolve *) params;

    double Ptot = theTsolve->rhoe_or_Ptot;
    double Sigma = theTsolve->Sigma;
    double Gamma = theTsolve->Gamma;
    double xi_r = theTsolve->xi_r;
    double Omega_eff = theTsolve->Omega_eff;

    double n = 1./(Gamma-1.);
    double cs = sqrt(Ptot/Sigma);
    double h = cs/Omega_eff;


    *y = Ptot - Sigma*T - xi_r* h *pow(T,4);

    *dy = - Sigma - 4.*xi_r*h*pow(T,3);
}

