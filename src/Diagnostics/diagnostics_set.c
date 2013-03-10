#define DIAGNOSTICS_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include "../Headers/Sim.h"
#include "../Headers/Cell.h"
#include "../Headers/Diagnostics.h"
#include "../Headers/TimeStep.h"
#include "../Headers/header.h"

void diagnostics_set(struct Diagnostics * theDiagnostics,struct Cell *** theCells,struct Sim * theSim,struct TimeStep * theTimeStep){
  int num_r_points = sim_N_r(theSim)-sim_Nghost_rmin(theSim)-sim_Nghost_rmax(theSim);
  int num_r_points_global = sim_N_r_global(theSim);

  double * VectorDiag_temp = malloc(sizeof(double) * num_r_points_global*theDiagnostics->NUM_DIAG);
  double * VectorDiag_reduce = malloc(sizeof(double) * num_r_points_global*theDiagnostics->NUM_DIAG);
  double * ScalarDiag_temp = malloc(sizeof(double)*theDiagnostics->NUM_DIAG);
  double * ScalarDiag_reduce = malloc(sizeof(double)*theDiagnostics->NUM_DIAG);

  double dtout = timestep_get_t(theTimeStep)-theDiagnostics->toutprev;
  
  int imin = sim_Nghost_rmin(theSim);
  int imax = sim_N_r(theSim)-sim_Nghost_rmax(theSim);
  int kmin = sim_Nghost_zmin(theSim);
  int kmax = sim_N_z(theSim)-sim_Nghost_zmax(theSim);
  int i,j,k,n;
  
  for (n=0;n<theDiagnostics->NUM_DIAG;++n){
    ScalarDiag_temp[n]=0.0;
    ScalarDiag_reduce[n]=0.0;
  }
  for (i=0;i<num_r_points_global;++i){
    for (n=0;n<theDiagnostics->NUM_DIAG;++n){
      VectorDiag_temp[i*theDiagnostics->NUM_DIAG+n]=0.0;
    }
  }

  for (k=kmin;k<kmax;++k){
    double zp = sim_z_faces(theSim,k);
    double zm = sim_z_faces(theSim,k-1);
    double dz = zp-zm;
    for (i=imin;i<imax;++i){
     for (j=0;j<sim_N_p(theSim,i);++j){
        // divide by number of phi cells to get phi average, mult by dz because we are doing a z integration;
        VectorDiag_temp[(sim_N_r_0(theSim)+i-imin)*theDiagnostics->NUM_DIAG+0] += (cell_prim(cell_single(theCells,i,j,k),RHO)/sim_N_p(theSim,i)*dz) ;
        VectorDiag_temp[(sim_N_r_0(theSim)+i-imin)*theDiagnostics->NUM_DIAG+1] += (cell_prim(cell_single(theCells,i,j,k),PPP)/sim_N_p(theSim,i)*dz) ;
        //etc
      }
    }
  }

  for (i=imin;i<imax;++i){
    double rp = sim_r_faces(theSim,i);
    double rm = sim_r_faces(theSim,i-1);
    for (n=0;n<theDiagnostics->NUM_DIAG;++n){
      // mult by delta r^2 because we are doing an r integration
      ScalarDiag_temp[n] += VectorDiag_temp[(sim_N_r_0(theSim)+i-imin)*theDiagnostics->NUM_DIAG+n] * (rp*rp-rm*rm); 
    }
  }

  MPI_Allreduce( ScalarDiag_temp,ScalarDiag_reduce, theDiagnostics->NUM_DIAG, MPI_DOUBLE, MPI_SUM, sim_comm);
  MPI_Allreduce( VectorDiag_temp,VectorDiag_reduce, num_r_points_global*theDiagnostics->NUM_DIAG, MPI_DOUBLE, MPI_SUM, sim_comm);
  double RMIN = sim_RMIN(theSim);
  double RMAX = sim_RMAX(theSim);
  double ZMIN = sim_ZMIN(theSim);
  double ZMAX = sim_ZMAX(theSim);

  for (n=0;n<theDiagnostics->NUM_DIAG;++n){
    ScalarDiag_reduce[n] /= ((ZMAX-ZMIN)*(RMAX*RMAX-RMIN*RMIN));
  }

  for (i=0;i<num_r_points_global;++i){
    for (n=0;n<theDiagnostics->NUM_DIAG;++n){
      VectorDiag_reduce[i*theDiagnostics->NUM_DIAG+n] /=(ZMAX-ZMIN);
    }
  }

  //We are doing time averaged diagnostics, so mult by delta t and add it
  //We will divide by the total delta next time we save to disk;
  for (i=0;i<num_r_points_global;++i){
    for (n=0;n<theDiagnostics->NUM_DIAG;++n){
      theDiagnostics->VectorDiag[i][n] += VectorDiag_reduce[i*theDiagnostics->NUM_DIAG+n] * dtout;
    }
  }
  for (n=0;n<theDiagnostics->NUM_DIAG;++n){
    theDiagnostics->ScalarDiag[n] += ScalarDiag_reduce[n] * dtout;
  }

  //update output time;
  theDiagnostics->toutprev = timestep_get_t(theTimeStep);

  free(ScalarDiag_temp);
  free(ScalarDiag_reduce);
  free(VectorDiag_temp);
  free(VectorDiag_reduce);
}

