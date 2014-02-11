#define SIM_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include "../Headers/Sim.h"
#include "../Headers/MPIsetup.h"
#include "../Headers/header.h"


void sim_alloc_arr(struct Sim * theSim, struct MPIsetup * theMPIsetup) {

  int N_noghost[2];
  N_noghost[R_DIR] = theSim->N_global[R_DIR]/mpisetup_dim_NumProcs(theMPIsetup,R_DIR);
  N_noghost[Z_DIR] = theSim->N_global[Z_DIR]/mpisetup_dim_NumProcs(theMPIsetup,Z_DIR);

  int Nghost_min[2];
  int Nghost_max[2];
  if (mpisetup_dim_MyProc(theMPIsetup,0) == 0 && sim_NoInnerBC(theSim) == 1){
    Nghost_min[R_DIR] = 1;
  }else{
    Nghost_min[R_DIR]  = theSim->ng;
  }
  Nghost_max[R_DIR] = theSim->ng;
  if (sim_N_global(theSim,Z_DIR)==1){
    Nghost_min[Z_DIR]=0;
    Nghost_max[Z_DIR]=0;
  }else{
    Nghost_min[Z_DIR] = theSim->ng;
    Nghost_max[Z_DIR] = theSim->ng;
  }

  int N_withghost[2];
  N_withghost[R_DIR] = N_noghost[R_DIR]+Nghost_min[R_DIR]+Nghost_max[R_DIR];
  N_withghost[Z_DIR] = N_noghost[Z_DIR]+Nghost_min[Z_DIR]+Nghost_max[Z_DIR];

  theSim->N_noghost[R_DIR] = N_noghost[R_DIR];
  theSim->N_noghost[Z_DIR] = N_noghost[Z_DIR];

  theSim->Nghost_min[R_DIR] = Nghost_min[R_DIR];
  theSim->Nghost_max[R_DIR] = Nghost_max[R_DIR];
  theSim->Nghost_min[Z_DIR] = Nghost_min[Z_DIR];
  theSim->Nghost_max[Z_DIR] = Nghost_max[Z_DIR];

  theSim->N_p = (int *) malloc(N_withghost[R_DIR]*sizeof(int));
  theSim->r_faces   = (double *) malloc( (N_withghost[R_DIR]+1)*sizeof(double) );
  theSim->z_faces   = (double *) malloc( (N_withghost[Z_DIR]+1)*sizeof(double) );
}


