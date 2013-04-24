#define SIM_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include "../Headers/Sim.h"
#include "../Headers/header.h"
#include <math.h>

double sim_W_A(struct Sim * theSim,double r){
  if (sim_W_A_TYPE(theSim)==A_NONE){
    //no w_analytic
    return(0.0);
  }else if (sim_W_A_TYPE(theSim)==A_OMEGA20){
    //rigid rotation with Omega = 20
    return(r*20); 
  } else if (sim_W_A_TYPE(theSim)==A_KEPLER){
    //keplerian rotation
    return(pow(r,-0.5));
   } else if (sim_W_A_TYPE(theSim)==A_WEIRD){
    //weird rotation
    return(pow(r,-0.5)*exp(-.1/r));
  } else{
    printf("ERROR\n");
    exit(0);
  }
}

double sim_OM_A_DERIV(struct Sim * theSim,double r){
  if (sim_W_A_TYPE(theSim)==A_NONE){
    //no w_analytic
    return(0.0);
  }else if (sim_W_A_TYPE(theSim)==A_OMEGA20){
    //rigid rotation with Omega = 20
    return(0.0); 
  } else if (sim_W_A_TYPE(theSim)==A_KEPLER){
    //keplerian rotation
    return(-1.5*pow(r,-2.5));
  } else if (sim_W_A_TYPE(theSim)==A_WEIRD){
    //weird rotation
    return((.1/r/r-1.5/r)*pow(r,-1.5)*exp(-.1/r));
  } else{
    printf("ERROR\n");
    exit(0);
  }
}

