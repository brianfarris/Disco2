#define GRID_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include "../Headers/Grid.h"
#include "../Headers/header.h"


int grid_N_p(struct Grid *theGrid,int i){
  return(theGrid->N_p[i]);
}
double grid_r_faces(struct Grid *theGrid,int i){
  return(theGrid->r_faces[i+1]);
}
double grid_z_faces(struct Grid *theGrid,int k){
  return(theGrid->z_faces[k+1]);
}
int grid_N_r(struct Grid *theGrid){
  return(theGrid->N_r);
}
int grid_N_z(struct Grid *theGrid){
  return(theGrid->N_z);
}
int grid_Ncells(struct Grid *theGrid){
  return(theGrid->Ncells);
}
int grid_Ncells_global(struct Grid *theGrid){
  return(theGrid->Ncells_global);
}
int grid_offset(struct Grid *theGrid){
  return(theGrid->offset);
}
int grid_Nghost_rmin(struct Grid *theGrid){
  return(theGrid->Nghost_rmin);
}
int grid_Nghost_rmax(struct Grid *theGrid){
  return(theGrid->Nghost_rmax);
}
int grid_Nghost_zmin(struct Grid *theGrid){
  return(theGrid->Nghost_zmin);
}
int grid_Nghost_zmax(struct Grid *theGrid){
  return(theGrid->Nghost_zmax);
}
int grid_ng(struct Grid *theGrid){
  return(theGrid->ng);
}
int grid_N_z_global(struct Grid *theGrid){
  return(theGrid->N_z_global);
}


