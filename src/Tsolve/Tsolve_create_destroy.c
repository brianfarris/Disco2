#define TSOLVE_PRIVATE_DEFS
#include <stdlib.h>
#include "../Headers/Tsolve.h"

struct Tsolve * tsolve_create(double rhoe_or_Ptot, double Sigma,
        double Gamma, double xi_r, double Omega_eff,double guess){
  struct Tsolve *theTsolve = malloc(sizeof(struct Tsolve));
  theTsolve->rhoe_or_Ptot = rhoe_or_Ptot;
  theTsolve->Sigma = Sigma;
  theTsolve->Gamma = Gamma;
  theTsolve->xi_r = xi_r;
  theTsolve->Omega_eff = Omega_eff;
  
  theTsolve->guess = guess;
  return(theTsolve);
}

void tsolve_destroy(struct Tsolve * theTsolve){
  free(theTsolve);
}

