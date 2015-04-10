#define EOS_PRIVATE_DEFS
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "../Headers/header.h"
#include "../Headers/EOS.h"
#include "../Headers/Metric.h"
#include "../Headers/Sim.h"

double eos_qdot_pair_itoh(double rho, double T);

double eos_cool_none(double *prim, double H, struct Sim *theSim)
{
    return 0.0;
}

double eos_cool_isotherm(double *prim, double H, struct Sim *theSim)
{
    double T0 = sim_CoolPar1(theSim);
    double tau = sim_CoolPar2(theSim);

    return (prim[PPP]/prim[RHO] - T0) / tau;
}

double eos_cool_bb_es(double *prim, double H, struct Sim *theSim)
{
    //TODO: Work in 2D and 3D?
    double rho, T;
    if(sim_Background(theSim) != GRDISC)
    {
        double P;
        rho = (prim[RHO]/H)*eos_rho_scale;
        P = (prim[PPP]/H)*eos_rho_scale;
        T = P / rho * eos_mp*eos_c*eos_c; //in ergs
    }
    else
    {
        rho = prim[RHO]*eos_rho_scale;
        T = prim[TTT] * eos_mp*eos_c*eos_c; //in ergs
    }
    double kappa = 0.2;
    double h = H*eos_r_scale;
    double q0 = sim_CoolPar1(theSim);
    double q = q0 * 8.0*eos_sb * T*T*T*T/(3.0*kappa*rho*h);
    double Q = q / (eos_c*eos_c*eos_c * eos_rho_scale);

    //if(sim_NUM_Q(theSim) > sim_NUM_C(theSim))
    //    Q *= prim[sim_NUM_C(theSim)];

    return Q;
}

double eos_cool_bb_ff(double *prim, double H, struct Sim *theSim)
{
    //TODO: Work in 2D and 3D?
    double rho, T;
    if(sim_Background(theSim) != GRDISC)
    {
        double P;
        rho = (prim[RHO]/H)*eos_rho_scale;
        P = (prim[PPP]/H)*eos_rho_scale;
        T = P / rho * eos_mp*eos_c*eos_c; //in ergs
    }
    else
    {
        rho = prim[RHO]*eos_rho_scale;
        T = prim[TTT] * eos_mp*eos_c*eos_c; //in ergs
    }
    double kappa = 6.6e22 * rho * pow(T/eos_k, -3.5);
    double h = H*eos_r_scale;
    double q0 = sim_CoolPar1(theSim);
    double q = q0 * 8.0*eos_sb * T*T*T*T/(3.0*kappa*rho*h);
    double Q = q / (eos_c*eos_c*eos_c * eos_rho_scale);

    return Q;
}

double eos_cool_neutrino_aprx(double *prim, double H, struct Sim *theSim)
{
    double rho10, t11;
    double x1 = sim_CoolPar1(theSim);
    double xnuc = 1.0;

    if(sim_EOSType(theSim) == EOS_PWF)
    {
        xnuc = eos_xnuc_pwf(prim[RHO]*eos_rho_scale, 
                            prim[TTT]*eos_mp*eos_c*eos_c/eos_k);
    }

    if(sim_Background(theSim) != GRDISC)
    {
        double P;
        rho10 = (prim[RHO]/H)*eos_rho_scale / (1.0e10);
        P = (prim[PPP]/H)*eos_rho_scale;
        t11 = P / prim[RHO] * eos_mp*eos_c*eos_c/(1.0e11*eos_k);
    }
    else
    {
        rho10 = prim[RHO] * eos_rho_scale / (1.0e10);
        t11 = prim[TTT] * eos_mp*eos_c*eos_c/(1.0e11*eos_k);
    }
    double h = H * eos_r_scale;
    double q = (5.0e33*pow(t11,9) + 9.0e33*xnuc*rho10*pow(t11,6)) * h;
    double Q = q / (eos_c*eos_c*eos_c * eos_rho_scale);

    return x1*Q;
}

double eos_cool_neutrino_itoh(double *prim, double H, struct Sim *theSim)
{
    double rho, tk;
    double x1 = sim_CoolPar1(theSim); //Pair coefficient
    double x2 = sim_CoolPar2(theSim); //URCA coefficient
    double xnuc;

    if(sim_Background(theSim) != GRDISC)
    {
        double P;
        rho = (prim[RHO]/H)*eos_rho_scale;
        P = (prim[PPP]/H)*eos_rho_scale;
        tk = P / prim[RHO] * eos_mp*eos_c*eos_c/eos_k;
    }
    else
    {
        rho = prim[RHO] * eos_rho_scale;
        tk = prim[TTT] * eos_mp*eos_c*eos_c/eos_k;
    }
    
    if(sim_EOSType(theSim) == EOS_PWF)
        xnuc = eos_xnuc_pwf(rho, tk);
    else
        xnuc = 1.0;

    double h = H * eos_r_scale;
    double qpair = eos_qdot_pair_itoh(rho, tk);
    double qurca =  9.0e33 * xnuc * (rho*1.0e-10) * pow(tk*1.0e-11,6) * h;
    double Q = (x1*qpair + x2*qurca) / (eos_c*eos_c*eos_c * eos_rho_scale);

    return Q;
}

double eos_qdot_pair_itoh(double rho, double T)
{
    // Pair neutrino loss rate from Itoh et. al. 1989 (equiv 1996)
    // Equations will be marked from this version.
    // Input: rho - density in g/cm^3, T - temperature in Kelvin
    // Output: energy loss rate in erg/(cm^3 s)

    double mu = 2.0; //nucleons per electron

    //Vector and Axial couplings (eq. 2 & 3)
    double cv = 0.5 + 2.0*eos_sin2thw;
    double cv1 = 1.0 - cv;
    double ca = 0.5;
    double ca1 = 1.0 - ca;

    //My own shorthand for the coupling for each channel.
    double alpha  = cv*cv + ca*ca + eos_n_nu*(cv1*cv1 + ca1*ca1);
    double alpha1 = cv*cv - ca*ca + eos_n_nu*(cv1*cv1 - ca1*ca1);

    //Dimensionless T & rho (eq. 8 & 7)
    double la = T / 5.9302e9;
    double xi = pow(rho/(mu*1.0e9), 1.0/3.0) / la;

    //Fitting coefficients (Tb. 3)
    double a0, a1, a2, b1, b2, b3, c;
    a0 = 6.002e19;
    a1 = 2.084e20;
    a2 = 1.872e21;
    if(T < 1.0e10)
    {
        b1 = 9.383e-1;
        b2 = -4.141e-1;
        b3 = 5.829e-2;
        c = 5.5924;
    }
    else
    {
        b1 = 1.2383;
        b2 = -0.8141;
        b3 = 0.0;
        c = 4.9924;
    }

    //Eq. 20
    double f = la*la*la*(a0 + a1*xi + a2*xi*xi) * exp(-c*xi)
                 / (la*la*la*xi*xi*xi + b1*la*la + b2*la + b3);
    //Eq. 21
    double g = 1.0 - 13.04*la*la + 133.5*la*la*la*la + 1543*la*la*la*la*la*la
                + 918.6*la*la*la*la*la*la*la*la;
    //Eq. 19
    double q = pow(1.0 + rho / (mu*(7.692e7*la*la*la+9.715e6*sqrt(la))), -0.3) 
                / (10.7480*la*la + 0.3967*sqrt(la) + 1.0050);

    //Eq. 18
    double Q = 0.5 * (alpha + alpha1*q) * g * f * exp(-2.0/la);

    return Q;
}
