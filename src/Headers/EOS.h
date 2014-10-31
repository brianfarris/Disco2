#ifndef EOS_H
#define EOS_H

struct Metric;
struct Sim;

#ifdef EOS_PRIVATE_DEFS
static double eos_cool_isotherm_T;
#endif

//Initialize
void eos_init(struct Sim *);

//Cooling
double (*eos_cool)(double *, double *, struct Sim *);
double eos_cool_none(double *x, double *prim, struct Sim *theSim);
double eos_cool_isotherm(double *x, double *prim, struct Sim *theSim);
double eos_cool_bb_es(double *x, double *prim, struct Sim *theSim);
double eos_cool_visc(double *x, double *prim, struct Sim *theSim);

#endif
