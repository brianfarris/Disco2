#define SIM_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include "../Headers/Sim.h"
#include "../Headers/MPIsetup.h"
#include "../Headers/header.h"


void sim_alloc_arr(struct Sim * theSim, struct MPIsetup * theMPIsetup) {

  int N_r_noghost = theSim->N_r_global/mpisetup_dim_NumProcs(theMPIsetup)[0];
  int N_z_noghost = theSim->N_z_global/mpisetup_dim_NumProcs(theMPIsetup)[1];

  printf("mpisetup_dim_MyProc(theMPIsetup,0): %d\n",mpisetup_dim_MyProc(theMPIsetup,0));
  int Nghost_rmin;
  if (mpisetup_dim_MyProc(theMPIsetup,0)==0){
    Nghost_rmin  = 1;
  }else{
    Nghost_rmin  = theSim->ng;
  }
  int Nghost_rmax = theSim->ng;
  int Nghost_zmin = theSim->ng;
  int Nghost_zmax = theSim->ng;


  int N_r_withghost = N_r_noghost+Nghost_rmin+Nghost_rmax;
  int N_z_withghost = N_z_noghost+Nghost_zmin+Nghost_zmax;

  theSim->N_r_noghost = N_r_noghost;
  theSim->N_z_noghost = N_z_noghost;

  theSim->Nghost_rmin = Nghost_rmin;
  theSim->Nghost_rmax = Nghost_rmax;
  theSim->Nghost_zmin = Nghost_zmin;
  theSim->Nghost_zmax = Nghost_zmax;

  theSim->N_p = (int *) malloc(N_r_withghost*sizeof(int));
  theSim->r_faces   = (double *) malloc( (N_r_withghost+1)*sizeof(double) );
  theSim->z_faces   = (double *) malloc( (N_z_withghost+1)*sizeof(double) );
}


