#define SIM_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include "../Headers/Sim.h"
#include "../Headers/header.h"
#include <math.h>

double sim_rOm_a(struct Sim * theSim,double r,double a){
  if (sim_W_A_TYPE(theSim)==A_FIXED){
    //no w_analytic
    return(0.0);
  } else if (sim_W_A_TYPE(theSim)==A_KEPLER){
    //keplerian rotation
    return(pow(r,-0.5));
  }else if (sim_W_A_TYPE(theSim)==A_OMEGA20){
    //rigid rotation with Omega = 20
    return(r*0); 
  }else if (sim_W_A_TYPE(theSim)==A_MILOS){
    //Useful for Milos/MacFadyen disks
    //return((1.-exp(-pow((r/0.5),w_a_milos_index)))/sqrt(r)); 
    return(r*(pow(r,w_a_milos_index-1.5)+pow(a,w_a_milos_index-1.5))/(pow(r,w_a_milos_index)+pow(a,w_a_milos_index)));
  } else if (sim_W_A_TYPE(theSim)==A_POWERLAW){
    //Useful for a torus
    return(pow(r,-1.0));
  } else{
    printf("ERROR\n");
    exit(0);
  }
}

double sim_rdrOm_a(struct Sim * theSim,double r,double a){
  if (sim_W_A_TYPE(theSim)==A_FIXED){
    //no w_analytic
    return(0.0);
  } else if (sim_W_A_TYPE(theSim)==A_KEPLER){
    //keplerian rotation
    return(-1.5*pow(r,-1.5));
  }else if (sim_W_A_TYPE(theSim)==A_OMEGA20){
    //rigid rotation with Omega = 20
    return(0.0); 
  }else if (sim_W_A_TYPE(theSim)==A_MILOS){
    //Useful for Milos/MacFadyen disks
    //return(1.5*pow(r,-1.5)*(-1.+(1.+(w_a_milos_index/1.5)*pow((r/0.5),w_a_milos_index))*exp(-pow((r/0.5),w_a_milos_index)))); 
    int n = w_a_milos_index;
    double Omega_a = (pow(r,n-1.5) + pow(a,n-1.5))/(pow(r,n) + pow(a,n));
    double arn = pow(a/r,n);
    return((-1.5 + n*arn/(1.+arn))*Omega_a + (1.5-n)*pow(a/r,n-1.5)/(1.+arn)*pow(r,-1.5));
  } else if (sim_W_A_TYPE(theSim)==A_POWERLAW){
    //Useful for a torus
    return(-2.*pow(r,-2.0));
  } else{
    printf("ERROR\n");
    exit(0);
  }
}

double sim_dtOm_a(struct Sim * theSim,double r,double a){
  if (sim_W_A_TYPE(theSim)==A_MILOS){
    double a0 = 1.0;
    double tau = sim_OrbShrinkTscale(theSim);
    double dta_o_a = -.25 * pow(a/a0,-4) / tau;
    int n = w_a_milos_index;
    double Omega_a = (pow(r,n-1.5) + pow(a,n-1.5))/(pow(r,n) + pow(a,n));
    double arn = pow(a/r,n);
    if (a>0.0){
      return(dta_o_a * ( Omega_a*(-1.5 + n/(1. + arn) ) + (1.5 - n)*pow(r,-1.5)/(1. + arn)));
    } else{
      return(0.0);
    }
  }else{
    return(0.0);
  }
}
