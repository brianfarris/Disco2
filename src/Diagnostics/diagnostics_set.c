#define DIAGNOSTICS_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include "../Headers/Grid.h"
#include "../Headers/Cell.h"
#include "../Headers/Diagnostics.h"
#include "../Headers/header.h"

void diagnostics_set(struct Diagnostics * theDiagnostics, struct Cell *** theCells, struct Grid * theGrid){
  int NUM_DIAG = 2;
  int num_r_points = grid_N_r(theGrid)-grid_Nghost_rmin(theGrid)-grid_Nghost_rmax(theGrid);
  /*
  double * VectorDiag_temp = malloc(sizeof(double)*NUM_DIAG*num_r_points);
  double ScalarDiag_temp[NUM_DIAG];
  double ScalarDiag_reduce[NUM_DIAG];
  */

  
  double * VectorDiag_temp = malloc(sizeof(double) * num_r_points*NUM_DIAG);
  int i;
  double * ScalarDiag_temp = malloc(sizeof(double)*NUM_DIAG);
  double * ScalarDiag_reduce = malloc(sizeof(double)*NUM_DIAG);

  
  int imin = grid_Nghost_rmin(theGrid);
  int imax = grid_N_r(theGrid)-grid_Nghost_rmax(theGrid);
  int kmin = grid_Nghost_zmin(theGrid);
  int kmax = grid_N_z(theGrid)-grid_Nghost_zmax(theGrid);
  int j,k,n;
  
  for (n=0;n<NUM_DIAG;++n){
    ScalarDiag_temp[n]=0.0;
    ScalarDiag_reduce[n]=0.0;
  }
  for (i=imin;i<imax;++i){
    for (n=0;n<NUM_DIAG;++n){
      VectorDiag_temp[(i-imin)*NUM_DIAG+n]=0.0;
    }
  }

  for (k=kmin;k<kmax;++k){
    double zp = grid_z_faces(theGrid,k);
    double zm = grid_z_faces(theGrid,k-1);
    double dz = zp-zm;
    for (i=imin;i<imax;++i){
      double rp = grid_r_faces(theGrid,i);
      double rm = grid_r_faces(theGrid,i-1);
      for (j=0;j<grid_N_p(theGrid,i);++j){
        // divide by number of phi cells to get phi average, mult by dz because we are doing a z integration;
        VectorDiag_temp[(i-imin)*NUM_DIAG+0] += (cell_prim(cell_single(theCells,i,j,k),RHO)/grid_N_p(theGrid,i)*dz) ;
        VectorDiag_temp[(i-imin)*NUM_DIAG+1] += (cell_prim(cell_single(theCells,i,j,k),PPP)/grid_N_p(theGrid,i)*dz) ;
        //etc
      }
      for (n=0;n<NUM_DIAG;++n){
        // mult by delta r^2 because we are doing an r integration
        ScalarDiag_temp[n] += VectorDiag_temp[(i-imin)*NUM_DIAG+n] * (rp*rp-rm*rm); 
      }
    }
  }

  MPI_Allreduce( ScalarDiag_temp,ScalarDiag_reduce, NUM_DIAG, MPI_DOUBLE, MPI_SUM, grid_comm);
  double RMIN = grid_RMIN(theGrid);
  double RMAX = grid_RMAX(theGrid);
  double ZMIN = grid_ZMIN(theGrid);
  double ZMAX = grid_ZMAX(theGrid);

  for (n=0;n<NUM_DIAG;++n){
    ScalarDiag_reduce[n] /= ((ZMAX-ZMIN)*(RMAX*RMAX-RMIN*RMIN));
  }

  free(ScalarDiag_temp);
  free(ScalarDiag_reduce);
  free(VectorDiag_temp);
}

