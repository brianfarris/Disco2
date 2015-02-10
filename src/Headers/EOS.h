#ifndef EOS_H
#define EOS_H

struct Metric;
struct Sim;

#ifdef EOS_PRIVATE_DEFS

//All units in c.g.s.  Temperature measured in erg.
const static double eos_c = 2.99792458e10;      //Speed of Light
const static double eos_G = 6.6738e-8;          //Gravitational Constant
const static double eos_h = 6.62606957e-27;     //Planck Constant
const static double eos_k = 1.3806488e-16;      //Boltzmann (erg/K)
const static double eos_sb = 1.56055371e59;     //Stefan-Boltzmann
const static double eos_mp = 1.672621777e-24;   //Proton Mass
const static double eos_rho_scale = 1.0;
const static double eos_r_scale = 1.4766250385e5;

const static double eos_n_nu = 3.0;             // # of Neutrino species
const static double eos_sin2thw = 0.2319;       // Sin^2(theta_weinberg)
const static double eos_alpha_bind = 4.53346658e-5; // alpha binding energy
#endif

//Initialize
void eos_init(struct Sim *);

//Equations of State
double (*eos_cs2)(double *, struct Sim *);
double (*eos_rho)(double *, struct Sim *);
double (*eos_eps)(double *, struct Sim *);
double (*eos_ppp)(double *, struct Sim *);
double (*eos_dpppdrho)(double *, struct Sim *);
double (*eos_dpppdttt)(double *, struct Sim *);
double (*eos_depsdrho)(double *, struct Sim *);
double (*eos_depsdttt)(double *, struct Sim *);

double eos_cs2_gammalaw(double *, struct Sim *);
double eos_rho_gammalaw(double *, struct Sim *);
double eos_eps_gammalaw(double *, struct Sim *);
double eos_ppp_gammalaw(double *, struct Sim *);
double eos_dpppdrho_gammalaw(double *, struct Sim *);
double eos_dpppdttt_gammalaw(double *, struct Sim *);
double eos_depsdrho_gammalaw(double *, struct Sim *);
double eos_depsdttt_gammalaw(double *, struct Sim *);

double eos_cs2_gasrad(double *, struct Sim *);
double eos_rho_gasrad(double *, struct Sim *);
double eos_eps_gasrad(double *, struct Sim *);
double eos_ppp_gasrad(double *, struct Sim *);
double eos_dpppdrho_gasrad(double *, struct Sim *);
double eos_dpppdttt_gasrad(double *, struct Sim *);
double eos_depsdrho_gasrad(double *, struct Sim *);
double eos_depsdttt_gasrad(double *, struct Sim *);

double eos_cs2_gasraddeg(double *, struct Sim *);
double eos_rho_gasraddeg(double *, struct Sim *);
double eos_eps_gasraddeg(double *, struct Sim *);
double eos_ppp_gasraddeg(double *, struct Sim *);
double eos_dpppdrho_gasraddeg(double *, struct Sim *);
double eos_dpppdttt_gasraddeg(double *, struct Sim *);
double eos_depsdrho_gasraddeg(double *, struct Sim *);
double eos_depsdttt_gasraddeg(double *, struct Sim *);

double eos_cs2_pwf(double *, struct Sim *);
double eos_rho_pwf(double *, struct Sim *);
double eos_eps_pwf(double *, struct Sim *);
double eos_ppp_pwf(double *, struct Sim *);
double eos_dpppdrho_pwf(double *, struct Sim *);
double eos_dpppdttt_pwf(double *, struct Sim *);
double eos_depsdrho_pwf(double *, struct Sim *);
double eos_depsdttt_pwf(double *, struct Sim *);
double eos_xnuc_pwf(double, double);
double eos_dxnudrho_pwf(double, double, double);
double eos_dxnudttt_pwf(double, double, double);

//Cooling
double (*eos_cool)(double *, double, struct Sim *);
double eos_cool_none(double *, double, struct Sim *);
double eos_cool_isotherm(double *, double, struct Sim *);
double eos_cool_bb_es(double *, double, struct Sim *);
double eos_cool_bb_ff(double *, double, struct Sim *);
double eos_cool_neutrino_aprx(double *, double, struct Sim *);
double eos_cool_neutrino_itoh(double *, double, struct Sim *);

#endif
