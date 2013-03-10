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
    for( i=0 ; i<sim_N_r(theSim) ; ++i ){
      theSim->N_p[i]=theSim->NP_CONST;
    }
  }else{
    for( i=0 ; i<sim_N_r(theSim) ; ++i ){
      double r = sim_r_faces(theSim,i);
      double dr = sim_r_faces(theSim,i)-sim_r_faces(theSim,i-1);
      theSim->N_p[i] = (int)( 2.*M_PI*( 1. + (r/dr-1.)/theSim->aspect ) ) ;
    }
  }

}

void sim_set_rz(struct Sim * theSim,struct MPIsetup * theMPIsetup){

  int N_r_0 = theSim->N_r_noghost*mpisetup_dim_MyProc(theMPIsetup,0);
  int N_z_0 = theSim->N_z_noghost*mpisetup_dim_MyProc(theMPIsetup,1);

  theSim->N_r_0 = N_r_0;
  theSim->N_z_0 = N_z_0;

  int i,k;
  for(i = 0; i < sim_N_r(theSim)+1; i++){
    int ig = i-theSim->Nghost_rmin+N_r_0;
    double delta = (theSim->RMAX-theSim->RMIN)/(double)theSim->N_r_global;
    theSim->r_faces[i] = theSim->RMIN+(double)ig*delta;
  }
  for(k = 0; k < sim_N_z(theSim)+1; k++){
    int kg = k-theSim->Nghost_zmin+N_z_0;
    double delta = (theSim->ZMAX-theSim->ZMIN)/(double)theSim->N_z_global;
    theSim->z_faces[k] = theSim->ZMIN+(double)kg*delta;
  } 
}

void sim_set_misc(struct Sim *theSim,struct MPIsetup * theMPIsetup) {
  int i,j,k,q;

  // Stuff that involves counting cells  
  int Ncells=0;
  int Ncells_global;

  for (i=theSim->Nghost_rmin;i<theSim->N_r_noghost+theSim->Nghost_rmin;++i){
    Ncells += theSim->N_p[i]*theSim->N_z_noghost;
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




