#define GRID_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include "../Headers/Grid.h"
#include "../Headers/MPIsetup.h"
#include "../Headers/header.h"

void grid_set_N_p(struct Grid * theGrid){
  int i;
  for(i = 0; i < grid_N_r(theGrid); i++){
    theGrid->N_p[i] = 1024;
  }
}

void grid_set_rz(struct Grid * theGrid,struct MPIsetup * theMPIsetup){
 
  int N_r_0 = theGrid->N_r_noghost*mpisetup_dim_MyProc(theMPIsetup,0);
  int N_z_0 = theGrid->N_z_noghost*mpisetup_dim_MyProc(theMPIsetup,1);

  int i,k;
  for(i = 0; i < grid_N_r(theGrid)+1; i++){
    int ig = i-theGrid->Nghost_rmin+N_r_0;
    double delta = (theGrid->RMAX-theGrid->RMIN)/(double)theGrid->N_r_global;
    theGrid->r_faces[i] = theGrid->RMIN+(double)ig*delta;
  }
  for(k = 0; k < grid_N_z(theGrid)+1; k++){
    int kg = k-theGrid->Nghost_zmin+N_z_0;
    double delta = (theGrid->ZMAX-theGrid->ZMIN)/(double)theGrid->N_z_global;
    theGrid->z_faces[k] = theGrid->ZMIN+(double)kg*delta;
  } 
}

void grid_set_Ncells_and_offset(struct Grid *theGrid,struct MPIsetup * theMPIsetup) {
  int i,j,k,q;

  int Ncells=0;
  int Ncells_global;

  for (i=0;i<theGrid->N_r_noghost;++i){
    Ncells += theGrid->N_p[i]*theGrid->N_z_noghost;
  }

  int *Ncells_arr = malloc(sizeof(int) * mpisetup_NumProcs(theMPIsetup));
  MPI_Allgather(&Ncells, 1, MPI_INT, Ncells_arr, 1, MPI_INT, MPI_COMM_WORLD);

  int offset=0;
  for (i=0;i<mpisetup_MyProc(theMPIsetup);i++){
    offset += Ncells_arr[i];
  }
  free(Ncells_arr);

  MPI_Allreduce( &Ncells , &Ncells_global , 1 , MPI_INT , MPI_SUM , grid_comm );

  theGrid->Ncells = Ncells;
  theGrid->Ncells_global = Ncells_global;
  theGrid->offset = offset;
}




