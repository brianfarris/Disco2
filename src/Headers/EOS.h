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
#endif

//Initialize
void eos_init(struct Sim *);

//Equations of State
double (*eos_cs2_prim)(double *, double, struct Sim *);
double (*eos_eps_prim)(double *, double, struct Sim *);
double (*eos_temp_prim)(double *, double, struct Sim *);

double eos_cs2_prim_ideal_newt(double *prim, double H, struct Sim *theSim);
double eos_eps_prim_ideal_newt(double *prim, double H, struct Sim *theSim);
double eos_temp_prim_ideal_newt(double *prim, double H, struct Sim *theSim);

double eos_cs2_prim_ideal_gr(double *prim, double H, struct Sim *theSim);
double eos_eps_prim_ideal_gr(double *prim, double H, struct Sim *theSim);
double eos_temp_prim_ideal_gr(double *prim, double H, struct Sim *theSim);

//Cooling
double (*eos_cool)(double *, double, struct Sim *);
double eos_cool_none(double *prim, double, struct Sim *theSim);
double eos_cool_isotherm(double *prim, double, struct Sim *theSim);
double eos_cool_bb_es(double *prim, double, struct Sim *theSim);
double eos_cool_bb_ff(double *prim, double, struct Sim *theSim);
double eos_cool_neutrino(double *prim, double, struct Sim *theSim);

#endif
