#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Sim.h"
#include "../Headers/header.h"

double cell_W_A(struct Sim * theSim,double r){
  if (sim_W_A_TYPE(theSim)==A_NONE){
    //no w_analytic
    return(0.0);
  }else if (sim_W_A_TYPE(theSim)==A_OMEGA10){
    //rigid rotation with Omega = 10
    return(r*10); 
  } else if (sim_W_A_TYPE(theSim)==A_KEPLER){
    //keplerian rotation
    return(pow(r,-1.5));
  }
}

double cell_OM_A_DERIV(struct Sim * theSim,double r){
  if (sim_W_A_TYPE(theSim)==A_NONE){
    //no w_analytic
    return(0.0);
  }else if (sim_W_A_TYPE(theSim)==A_OMEGA10){
    //rigid rotation with Omega = 10
    return(0.0); 
  } else if (sim_W_A_TYPE(theSim)==A_KEPLER){
    //keplerian rotation
    return(-1.5*pow(r,-2.5));
  }
}


