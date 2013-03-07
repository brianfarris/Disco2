#define GRID_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include "../Headers/Grid.h"
#include "../Headers/MPIsetup.h"
#include "../Headers/header.h"


void grid_alloc_arr(struct Grid * theGrid, struct MPIsetup * theMPIsetup) {

  int N_r_noghost = theGrid->N_r_global/mpisetup_dim_NumProcs(theMPIsetup)[0];
  int N_z_noghost = theGrid->N_z_global/mpisetup_dim_NumProcs(theMPIsetup)[1];

  printf("mpisetup_dim_MyProc(theMPIsetup,0): %d\n",mpisetup_dim_MyProc(theMPIsetup,0));
  int Nghost_rmin;
  if (mpisetup_dim_MyProc(theMPIsetup,0)==0){
    Nghost_rmin  = 1;
  }else{
    Nghost_rmin  = theGrid->ng;
  }
  int Nghost_rmax = theGrid->ng;
  int Nghost_zmin = theGrid->ng;
  int Nghost_zmax = theGrid->ng;


  int N_r_withghost = N_r_noghost+Nghost_rmin+Nghost_rmax;
  int N_z_withghost = N_z_noghost+Nghost_zmin+Nghost_zmax;

  theGrid->N_r_noghost = N_r_noghost;
  theGrid->N_z_noghost = N_z_noghost;

  theGrid->Nghost_rmin = Nghost_rmin;
  theGrid->Nghost_rmax = Nghost_rmax;
  theGrid->Nghost_zmin = Nghost_zmin;
  theGrid->Nghost_zmax = Nghost_zmax;

  theGrid->N_p = (int *) malloc(N_r_withghost*sizeof(int));
  theGrid->r_faces   = (double *) malloc( (N_r_withghost+1)*sizeof(double) );
  theGrid->z_faces   = (double *) malloc( (N_z_withghost+1)*sizeof(double) );
}


