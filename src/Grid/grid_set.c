#define GRID_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include "../Headers/grid.h"
#include "../Headers/header.h"

void grid_set_N_p(struct Grid * theGrid){
  int i;
  for(i = 0; i < theGrid->N_r+theGrid->Nghost_rmin+theGrid->Nghost_rmax; i++){
    theGrid->N_p[i] = 16;
  }
}

void grid_set_rz(struct Grid * theGrid){
  int N_r_0 = theGrid->N_r*dim_MyProc[0];
  int N_z_0 = theGrid->N_z*dim_MyProc[1];

  int i;
  for(i = 0; i < theGrid->N_r+theGrid->Nghost_rmin+theGrid->Nghost_rmax+1; i++){
    int ig = i-theGrid->Nghost_rmin+N_r_0;
    double delta = (RMAX-RMIN)/(double)N_r_global;
    theGrid->r_faces[i] = RMIN+(double)ig*delta;
  }
  for(i = 0; i < theGrid->N_z+theGrid->Nghost_zmin+theGrid->Nghost_zmax+1; i++){
    int ig = i-theGrid->Nghost_zmin+N_z_0;
    double delta = (ZMAX-ZMIN)/(double)N_z_global;
    theGrid->z_faces[i] = ZMIN+(double)ig*delta;
  } 
}

void grid_set_Ncells_and_offset(struct Grid *theGrid) {
  int i,j,k,q;

  int Ncells=0;
  int Ncells_global;

  for (i=0;i<theGrid->N_r;++i){
    Ncells += theGrid->N_p[i]*theGrid->N_z;
  }

  int *Ncells_arr = malloc(sizeof(int) * NumProcs);
  MPI_Allgather(&Ncells, 1, MPI_INT, Ncells_arr, 1, MPI_INT, MPI_COMM_WORLD);

  int offset=0;
  for (i=0;i<MyProc;i++){
    offset += Ncells_arr[i];
  }
  free(Ncells_arr);

  MPI_Allreduce( &Ncells , &Ncells_global , 1 , MPI_INT , MPI_SUM , grid_comm );

  theGrid->Ncells = Ncells;
  theGrid->Ncells_global = Ncells_global;
  theGrid->offset = offset;
}



