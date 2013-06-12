#define SIM_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Sim.h"
#include "../Headers/Cell.h"
#include "../Headers/MPIsetup.h"
#include "../Headers/header.h"

void sim_set_N_p(struct Sim * theSim){
  int i;
  if (theSim->NP_CONST>0){
    for( i=0 ; i<sim_N(theSim,R_DIR) ; ++i ){
      theSim->N_p[i]=theSim->NP_CONST;
    }
  }else{
    for( i=0 ; i<sim_N(theSim,R_DIR) ; ++i ){
      double r = sim_FacePos(theSim,i,R_DIR);
      double dr = sim_FacePos(theSim,i,R_DIR)-sim_FacePos(theSim,i-1,R_DIR);
      theSim->N_p[i] = (int)( 2.*M_PI*( 1. + (r/dr-1.)/theSim->aspect ) ) ;
      if (theSim->NPCAP>0){
        if (theSim->N_p[i]>theSim->NPCAP) theSim->N_p[i] = theSim->NPCAP;
      }
    }
  }
}

//used by root finder
double r_func( double r, double r2, double fac,double sigma, double r0){
  double F = r+fac*sigma*sqrt(M_PI/4.0)*(erf((r-r0)/sigma)+1.0)-r2;
  return(F);
}

//find the sign
int sgn(double x) {
  if (x > 0) return 1;
  if (x < 0) return -1;
  return 0;
}

//root finder
double get_r(double r2,double RMIN,double RMAX,double r0,double sigma,double fac){
  double r;
  int n=1;
  int NMAX = 1000;
  double a=RMIN;
  double b=RMAX;
  double TOL=1.0e-12;
  while (n <= NMAX) { 
    double c = (a + b)/2.;
    if ((r_func(c,r2,fac,sigma,r0) == 0.0)||((b-a)/2.<TOL)){
      r=c;
      break;
    }
    n=n+1;
    if (sgn(r_func(c,r2,fac,sigma,r0)) == sgn(r_func(a,r2,fac,sigma,r0))){
      a = c;
    }else {
      b = c;
    }
  }
  return(r);
}

void sim_set_rz(struct Sim * theSim,struct MPIsetup * theMPIsetup){

  int N0[2];
  N0[R_DIR] = theSim->N_noghost[R_DIR]*mpisetup_dim_MyProc(theMPIsetup,0);
  N0[Z_DIR] = theSim->N_noghost[Z_DIR]*mpisetup_dim_MyProc(theMPIsetup,1);

  theSim->N0[R_DIR] = N0[R_DIR];
  theSim->N0[Z_DIR] = N0[Z_DIR];

  double RMIN = theSim->MIN[R_DIR];
  double ZMIN = theSim->MIN[Z_DIR];
  double RMAX=theSim->MAX[R_DIR];
  double ZMAX=theSim->MAX[Z_DIR];

  //****************************************************************************************************************************
  // Here, we are setting up the distribution of radial points. We would like to have two freedoms:
  // 1) Space them logarithmically so that there is higher resolution near r=0 and lower resolution far away.
  //    We would like to control the scale over which the spacing gets logarithmic, so that for r<<Rscale, the spacing is uniform.
  // 2) Have the ability to increase the resolution around some other radius r=r0. 
  //    This could, for example, correspond to the radius r0 at which a binary is located. 
  //    We would like the resolution to improve at r0 by a factor of 1+fac, and go back to normal over a scale of sigma.
  //****************************************************************************************************************************

  // Rscale is the scale at which the points begin to get logarithmically spaced
  double Rscale = theSim->RLogScale; 
  //sigma is the approx size of the hi-res region
  double sigma = theSim->HiResSigma;
  // r0 is the radius at which we center a hi-res region
  double r0 = theSim->HiResR0;
  // how much the res increases. res changes by factor of 1+fac at r0.
  double fac = theSim->HiResFac;

  // r1 is a different coordinate in which the cells are evenly spaced. Here we find the max/min of r1.
  double r2_max = RMAX + fac*sigma*sqrt(M_PI/4.0)*(erf((RMAX-r0)/sigma)+1.0);
  double r2_min = RMIN + fac*sigma*sqrt(M_PI/4.0)*(erf((RMIN-r0)/sigma)+1.0);
  double r1_max = Rscale*log(1.0+r2_max/Rscale);
  double r1_min = Rscale*log(1.0+r2_min/Rscale);

  int i,k;
  for(i = 0; i < sim_N(theSim,R_DIR)+1; i++){
    int ig;
    ig = i-theSim->Nghost_min[R_DIR]+N0[R_DIR];
    double delta = (r1_max-r1_min)/(double)theSim->N_global[R_DIR];
    double r1 = r1_min+(double)ig*delta;
    double r2 = Rscale*(exp(r1/Rscale)-1.0);

    theSim->r_faces[i] = get_r(r2,0.,2.*RMAX,r0,sigma,fac);

  }

  // if we have no inner BC, then the inner zones should extend to the origin
  if ((theSim->NoInnerBC == 1)&&(mpisetup_check_rin_bndry(theMPIsetup))){
    theSim->r_faces[0] = 0.0;
  }

  // Zscale is the scale at which the points begin to get logarithmically spaced
  double Zscale = theSim->ZLogScale; 

  double z1_max = sgn(ZMAX)*Zscale*log(1.0+fabs(ZMAX)/Zscale);
  double z1_min = sgn(ZMIN)*Zscale*log(1.0+fabs(ZMIN)/Zscale);

  for(k = 0; k < sim_N(theSim,Z_DIR)+1; k++){
    int kg = k-theSim->Nghost_min[Z_DIR]+N0[Z_DIR];
    double delta = (z1_max-z1_min)/(double)theSim->N_global[Z_DIR];
    double z1 = z1_min+(double)kg*delta;
    theSim->z_faces[k] = sgn(z1)*Zscale*(exp(fabs(z1)/Zscale)-1.0);
  } 
}

void sim_set_misc(struct Sim *theSim,struct MPIsetup * theMPIsetup) {
  int i,j,k,q;

  // Stuff that involves counting cells  
  int Ncells=0;
  int Ncells_global;

  //for (i=theSim->Nghost_min[R_DIR];i<theSim->N_noghost[R_DIR]+theSim->Nghost_min[R_DIR];++i){
  for (i=0;i<sim_N(theSim,R_DIR);++i){
    Ncells += theSim->N_p[i]*sim_N(theSim,Z_DIR);
  }

  int *Ncells_arr = malloc(sizeof(int) * mpisetup_NumProcs(theMPIsetup));
  MPI_Allgather(&Ncells, 1, MPI_INT, Ncells_arr, 1, MPI_INT, MPI_COMM_WORLD);

  int offset=0;
  for (i=0;i<mpisetup_MyProc(theMPIsetup);i++){
    offset += Ncells_arr[i];
  }
  free(Ncells_arr);

  MPI_Allreduce( &Ncells , &Ncells_global , 1 , MPI_INT , MPI_SUM , sim_comm );

  theSim->Ncells = Ncells;
  theSim->Ncells_global = Ncells_global;
  theSim->offset = offset;

  // For now, we will always say that there are two masses, and we will set masses to 0 when we need to
  theSim->NumGravMass = 2;

  theSim->NUM_Q = theSim->NUM_C + theSim->NUM_N;

}




