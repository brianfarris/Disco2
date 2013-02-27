#define GRID_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include "../Headers/Grid.h"
#include "../Headers/MPIsetup.h"
#include "../Headers/header.h"


struct Grid *grid_create(struct MPIsetup * theMPIsetup) {
  struct Grid *theGrid = (struct Grid *) malloc(sizeof(struct Grid));

  /*
  theGrid->Restart = 0;
  theGrid->N_r_global = 32;
  theGrid->N_z_global = 32;

  theGrid->ng = 2;

  theGrid->NUM_Q = 9;

  theGrid->Z_PERIODIC = 1;
  theGrid->MOVE_CELLS=C_RIGID;
  theGrid->NumGravMass = 1;

  theGrid->GAMMALAW=1.66666666;
  theGrid->INCLUDE_VISCOSITY=0;
  theGrid->EXPLICIT_VISCOSITY=0.1;
  theGrid->DIVB_CH=0.06;
  theGrid->DIVB_L=0.1;
  theGrid->CFL=0.5;
  theGrid->PLM=1.0;
  theGrid->POWELL=1;
  theGrid->GRAV2D=1;
  theGrid->G_EPS=0.0;
  theGrid->PHI_ORDER=2.0;
  theGrid->RHO_FLOOR=0.00001;
  theGrid->CS_FLOOR=0.0001;
  theGrid->CS_CAP=1.0;
  theGrid->VEL_CAP=10.0;
  theGrid->T_MAX = 50.0;
  theGrid->NUM_CHECKPOINTS = 100;
  */
  grid_read_par_file(theGrid,theMPIsetup,"mri_flock.par");

  int N_r = theGrid->N_r_global/mpisetup_dim_NumProcs(theMPIsetup)[0];
  int N_z = theGrid->N_z_global/mpisetup_dim_NumProcs(theMPIsetup)[1];

  int Nghost_rmin;
  if (mpisetup_dim_MyProc(theMPIsetup)[0]==0){
    Nghost_rmin  = 1;
  }else{
    Nghost_rmin  = theGrid->ng;
  }
  int Nghost_rmax = theGrid->ng;
  int Nghost_zmin = theGrid->ng;
  int Nghost_zmax = theGrid->ng;

  int N_r_withghost = N_r+Nghost_rmin+Nghost_rmax;
  int N_z_withghost = N_z+Nghost_zmin+Nghost_zmax;


  theGrid->N_z = N_z;
  theGrid->N_r = N_r;

  theGrid->Nghost_rmin = Nghost_rmin;
  theGrid->Nghost_rmax = Nghost_rmax;
  theGrid->Nghost_zmin = Nghost_zmin;
  theGrid->Nghost_zmax = Nghost_zmax;

  theGrid->N_p = (int *) malloc(N_r_withghost*sizeof(int));
  theGrid->r_faces   = (double *) malloc( (N_r_withghost+1)*sizeof(double) );
  theGrid->z_faces   = (double *) malloc( (N_z_withghost+1)*sizeof(double) );
  return theGrid;
}

void grid_destroy(struct Grid *theGrid) {
  free(theGrid->r_faces);
  free(theGrid->z_faces);
  free(theGrid->N_p);
  free(theGrid);
}


