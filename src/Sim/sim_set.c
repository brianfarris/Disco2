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
    }
  }

}

void sim_set_rz(struct Sim * theSim,struct MPIsetup * theMPIsetup){

  int N0[2];
  N0[R_DIR] = theSim->N_noghost[R_DIR]*mpisetup_dim_MyProc(theMPIsetup,0);
  N0[Z_DIR] = theSim->N_noghost[Z_DIR]*mpisetup_dim_MyProc(theMPIsetup,1);

  theSim->N0[R_DIR] = N0[R_DIR];
  theSim->N0[Z_DIR] = N0[Z_DIR];

  int i,k;
  for(i = 0; i < sim_N(theSim,R_DIR)+1; i++){
    int ig = i-theSim->Nghost_min[R_DIR]+N0[R_DIR];
    double delta = (theSim->MAX[R_DIR]-theSim->MIN[R_DIR])/(double)theSim->N_global[R_DIR];
    theSim->r_faces[i] = theSim->MIN[R_DIR]+(double)ig*delta;
  }
  for(k = 0; k < sim_N(theSim,Z_DIR)+1; k++){
    int kg = k-theSim->Nghost_min[Z_DIR]+N0[Z_DIR];
    double delta = (theSim->MAX[Z_DIR]-theSim->MIN[Z_DIR])/(double)theSim->N_global[Z_DIR];
    theSim->z_faces[k] = theSim->MIN[Z_DIR]+(double)kg*delta;
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

}




