#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Grid.h"
#include "../Headers/Face.h"
#include "../Headers/GravMass.h"
#include "../Headers/MPIsetup.h"
#include "../Headers/header.h"

struct Cell ***cell_create(struct Grid *theGrid,struct MPIsetup * theMPIsetup){
  int N_r_withghost = grid_N_r(theGrid)+grid_Nghost_rmin(theGrid)+grid_Nghost_rmax(theGrid);
  int N_z_withghost = grid_N_z(theGrid)+grid_Nghost_zmin(theGrid)+grid_Nghost_zmax(theGrid);

  //count up the total number of cells
  int Ncells_withghost=0;
  int i,j,k;
  for (k=0; k<N_z_withghost; k++) {
    for (i=0; i<N_r_withghost;i++){
      for(j = 0; j < grid_N_p(theGrid,i); j++){
        ++Ncells_withghost;
      }
    }
  }

  //allocate memory contiguously
  struct Cell ***theCells = (struct Cell ***)malloc(sizeof(struct Cell **)*N_z_withghost);
  theCells[0] = (struct Cell **)malloc(sizeof(struct Cell *)*N_r_withghost*N_z_withghost);
  theCells[0][0] = (struct Cell *)malloc(sizeof(struct Cell)*Ncells_withghost);
  
  for (k=0; k<N_z_withghost; k++) {
    theCells[k] = theCells[0]+k*N_r_withghost;
    for (i=0; i<N_r_withghost;i++){
      theCells[k][i] = theCells[0][0]+grid_N_p(theGrid,i)*(k*N_r_withghost+i);
    }
  }

  //  initialize
  srand(666+mpisetup_MyProc(theMPIsetup));
  for(k = 0; k < N_z_withghost; k++){
    for(i = 0; i < N_r_withghost; i++){
      double tiph_0 = 2.*M_PI*(double)rand()/(double)RAND_MAX;
      double dphi = 2.*M_PI/(double)grid_N_p(theGrid,i);
      for(j = 0; j < grid_N_p(theGrid,i); j++){
        double tiph = tiph_0 + (double)j*dphi;
        theCells[k][i][j].prim[RHO] = 0.0;
        theCells[k][i][j].prim[PPP] = 0.0;
        theCells[k][i][j].prim[URR] = 0.0;
        theCells[k][i][j].prim[UPP] = 0.0;
        theCells[k][i][j].prim[UZZ] = 0.0;
        theCells[k][i][j].prim[BRR] = 0.0;
        theCells[k][i][j].prim[BPP] = 0.0;
        theCells[k][i][j].prim[BZZ] = 0.0;
        theCells[k][i][j].prim[PSI] = 0.0;
        theCells[k][i][j].wiph = 0.0;
        theCells[k][i][j].divB = 0.0;
        theCells[k][i][j].GradPsi[0] = 0.0;
        theCells[k][i][j].GradPsi[1] = 0.0;
        theCells[k][i][j].GradPsi[2] = 0.0;
        theCells[k][i][j].tiph = tiph;
        theCells[k][i][j].RKtiph = tiph;
        theCells[k][i][j].dphi = dphi;
        theCells[k][i][j].wiph = 0.0;
        theCells[k][i][j].divB = 0.0;
      }
    }
  }

  return theCells;
}

void cell_destroy(struct Cell ***theCells,struct Grid *theGrid){
  int N_r_withghost = grid_N_r(theGrid)+grid_Nghost_rmin(theGrid)+grid_Nghost_rmax(theGrid);
  int N_z_withghost = grid_N_z(theGrid)+grid_Nghost_zmin(theGrid)+grid_Nghost_zmax(theGrid);
  int i,k;
  free(theCells[0][0]);
  free(theCells[0]);
  free(theCells);
}

