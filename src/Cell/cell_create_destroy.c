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
  int NUM_Q = grid_NUM_Q(theGrid);
  //count up the total number of cells
  int Ncells=0;
  int i,j,k;
  for (k=0; k<grid_N_z(theGrid); k++) {
    for (i=0; i<grid_N_r(theGrid);i++){
      for(j = 0; j < grid_N_p(theGrid,i); j++){
        ++Ncells;
      }
    }
  }

  //I wanted to allocate memory contiguously. It's not quite there yet.
  struct Cell ***theCells = (struct Cell ***)malloc(sizeof(struct Cell **)*grid_N_z(theGrid));
  theCells[0] = (struct Cell **)malloc(sizeof(struct Cell *)*grid_N_r(theGrid)*grid_N_z(theGrid));
  theCells[0][0] = (struct Cell *)malloc(sizeof(struct Cell)*Ncells);

  int position = 0;
  for (k=0; k<grid_N_z(theGrid); k++) {
    theCells[k] = theCells[0]+k*grid_N_r(theGrid);
    for (i=0; i<grid_N_r(theGrid);i++){
      theCells[k][i] = theCells[0][0]+position;//grid_N_p(theGrid,i)*(k*grid_N_r(theGrid)+i);
      position += grid_N_p(theGrid,i);
    }
  }
  Ncells=0;
  for (k=0; k<grid_N_z(theGrid); k++) {
    for (i=0; i<grid_N_r(theGrid);i++){
      for(j = 0; j < grid_N_p(theGrid,i); j++){
        theCells[k][i][j].prim = malloc(NUM_Q*sizeof(double));
        theCells[k][i][j].cons = malloc(NUM_Q*sizeof(double));
        theCells[k][i][j].RKcons = malloc(NUM_Q*sizeof(double));
        theCells[k][i][j].grad = malloc(NUM_Q*sizeof(double));
        theCells[k][i][j].gradp = malloc(NUM_Q*sizeof(double));
        ++Ncells;
      }
    }
  }

  //  initialize
  srand(666+mpisetup_MyProc(theMPIsetup));
  for(k = 0; k < grid_N_z(theGrid); k++){
    for(i = 0; i < grid_N_r(theGrid); i++){
      double tiph_0 = 2.*M_PI*(double)rand()/(double)RAND_MAX;
      double dphi = 2.*M_PI/(double)grid_N_p(theGrid,i);
      for(j = 0; j < grid_N_p(theGrid,i); j++){
        double tiph = tiph_0 + (double)j*dphi;
        theCells[k][i][j].prim[RHO] = 0.0;
        theCells[k][i][j].prim[PPP] = 0.0;
        theCells[k][i][j].prim[URR] = 0.0;
        theCells[k][i][j].prim[UPP] = 0.0;
        theCells[k][i][j].prim[UZZ] = 0.0;
        theCells[k][i][j].cons[DDD] = 0.0;
        theCells[k][i][j].cons[TAU] = 0.0;
        theCells[k][i][j].cons[SRR] = 0.0;
        theCells[k][i][j].cons[LLL] = 0.0;
        theCells[k][i][j].cons[SZZ] = 0.0;
        theCells[k][i][j].RKcons[DDD] = 0.0;
        theCells[k][i][j].RKcons[TAU] = 0.0;
        theCells[k][i][j].RKcons[SRR] = 0.0;
        theCells[k][i][j].RKcons[LLL] = 0.0;
        theCells[k][i][j].RKcons[SZZ] = 0.0;
        if (grid_runtype(theGrid)==1){
          theCells[k][i][j].prim[BRR] = 0.0;
          theCells[k][i][j].prim[BPP] = 0.0;
          theCells[k][i][j].prim[BZZ] = 0.0;
          theCells[k][i][j].prim[PSI] = 0.0;
          theCells[k][i][j].cons[BRR] = 0.0;
          theCells[k][i][j].cons[BPP] = 0.0;
          theCells[k][i][j].cons[BZZ] = 0.0;
          theCells[k][i][j].cons[PSI] = 0.0;
          theCells[k][i][j].RKcons[BRR] = 0.0;
          theCells[k][i][j].RKcons[BPP] = 0.0;
          theCells[k][i][j].RKcons[BZZ] = 0.0;
          theCells[k][i][j].RKcons[PSI] = 0.0;
        }
        theCells[k][i][j].wiph = 0.0;
        theCells[k][i][j].divB = 0.0;
        theCells[k][i][j].GradPsi[0] = 0.0;
        theCells[k][i][j].GradPsi[1] = 0.0;
        theCells[k][i][j].GradPsi[2] = 0.0;
        theCells[k][i][j].tiph = tiph;
        theCells[k][i][j].RKtiph = tiph;
        theCells[k][i][j].dphi = dphi;
      }
    }
  }

  return theCells;
}

void cell_destroy(struct Cell ***theCells,struct Grid *theGrid){
  int i,j,k;
  for (k=0; k<grid_N_z(theGrid); k++) {
    for (i=0; i<grid_N_r(theGrid);i++){
      for(j = 0; j < grid_N_p(theGrid,i); j++){
        free(theCells[k][i][j].gradp);
        free(theCells[k][i][j].grad);
        free(theCells[k][i][j].RKcons);
        free(theCells[k][i][j].cons);
        free(theCells[k][i][j].prim);
      }
    }
  }

  free(theCells[0][0]);
  free(theCells[0]);
  free(theCells);
}

