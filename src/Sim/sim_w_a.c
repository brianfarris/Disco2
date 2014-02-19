#define SIM_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include "../Headers/Sim.h"
#include "../Headers/header.h"
#include <math.h>

double sim_W_A(struct Sim * theSim,double r){
  if (sim_W_A_TYPE(theSim)==A_FIXED){
    //no w_analytic
    return(0.0);
  } else if (sim_W_A_TYPE(theSim)==A_KEPLER){
    //keplerian rotation
    return(pow(r,-0.5));
  }else if (sim_W_A_TYPE(theSim)==A_OMEGA20){
    //rigid rotation with Omega = 20
    return(r*20); 
  }else if (sim_W_A_TYPE(theSim)==A_MILOS){
    //Useful for Milos/MacFadyen disks
    return((1.-exp(-pow((r/0.5),1.5)))/sqrt(r)); 
  } else if (sim_W_A_TYPE(theSim)==A_POWERLAW){
    //Useful for a torus
    return(pow(r,-1.0));
  } else{
    printf("ERROR\n");
    exit(0);
  }
}

double sim_OM_A_DERIV(struct Sim * theSim,double r){
  if (sim_W_A_TYPE(theSim)==A_FIXED){
    //no w_analytic
    return(0.0);
  } else if (sim_W_A_TYPE(theSim)==A_KEPLER){
    //keplerian rotation
    return(-1.5*pow(r,-2.5));
  }else if (sim_W_A_TYPE(theSim)==A_OMEGA20){
    //rigid rotation with Omega = 20
    return(0.0); 
  }else if (sim_W_A_TYPE(theSim)==A_MILOS){
    //Useful for Milos/MacFadyen disks
    return(1.5*pow(r,-2.5)*(-1.+(1.+pow((r/0.5),1.5))*exp(-pow((r/0.5),1.5)))); 
  } else if (sim_W_A_TYPE(theSim)==A_POWERLAW){
    //Useful for a torus
    return(-2.*pow(r,-3.0));
  } else{
    printf("ERROR\n");
    exit(0);
  }
}

