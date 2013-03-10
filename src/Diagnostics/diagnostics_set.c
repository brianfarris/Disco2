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

  int NUM_SCAL = theDiagnostics->NUM_DIAG;
  int NUM_VEC = theDiagnostics->NUM_DIAG+1;

  double * VectorDiag_temp = malloc(sizeof(double) * num_r_points_global*NUM_VEC);
  double * VectorDiag_reduce = malloc(sizeof(double) * num_r_points_global*NUM_VEC);
  double * ScalarDiag_temp = malloc(sizeof(double)*NUM_SCAL);
  double * ScalarDiag_reduce = malloc(sizeof(double)*NUM_SCAL);

  double dtout = timestep_get_t(theTimeStep)-theDiagnostics->toutprev;
  
  int imin = sim_Nghost_rmin(theSim);
  int imax = sim_N_r(theSim)-sim_Nghost_rmax(theSim);
  int kmin = sim_Nghost_zmin(theSim);
  int kmax = sim_N_z(theSim)-sim_Nghost_zmax(theSim);
  int i,j,k,n;
  
  for (n=0;n<NUM_SCAL;++n){
    ScalarDiag_temp[n]=0.0;
    ScalarDiag_reduce[n]=0.0;
  }
  for (i=0;i<num_r_points_global;++i){
    for (n=0;n<NUM_VEC;++n){
      VectorDiag_temp[i*NUM_VEC+n]=0.0;
    }
  }

  for (k=kmin;k<kmax;++k){
    double zp = sim_z_faces(theSim,k);
    double zm = sim_z_faces(theSim,k-1);
    double dz = zp-zm;
    for (i=imin;i<imax;++i){
      double rp = sim_r_faces(theSim,i);
      double rm = sim_r_faces(theSim,i-1);
      double r = 0.5*(rm+rp);
      for (j=0;j<sim_N_p(theSim,i);++j){
        // divide by number of phi cells to get phi average, mult by dz because we are doing a z integration;
        VectorDiag_temp[(sim_N_r_0(theSim)+i-imin)*NUM_VEC+0] += r/sim_N_p(theSim,i)*dz ;
        VectorDiag_temp[(sim_N_r_0(theSim)+i-imin)*NUM_VEC+1] += (cell_prim(cell_single(theCells,i,j,k),RHO)/sim_N_p(theSim,i)*dz) ;
        VectorDiag_temp[(sim_N_r_0(theSim)+i-imin)*NUM_VEC+2] += (cell_prim(cell_single(theCells,i,j,k),PPP)/sim_N_p(theSim,i)*dz) ;
        //etc
      }
    }
  }

  for (i=imin;i<imax;++i){
    double rp = sim_r_faces(theSim,i);
    double rm = sim_r_faces(theSim,i-1);
    for (n=0;n<NUM_SCAL;++n){
      // mult by delta r^2 because we are doing an r integration
      ScalarDiag_temp[n] += VectorDiag_temp[(sim_N_r_0(theSim)+i-imin)*NUM_VEC+n+1] * (rp*rp-rm*rm); 
    }
  }

  MPI_Allreduce( ScalarDiag_temp,ScalarDiag_reduce , NUM_SCAL, MPI_DOUBLE, MPI_SUM, sim_comm);
  MPI_Allreduce( VectorDiag_temp,VectorDiag_reduce , num_r_points_global*NUM_VEC, MPI_DOUBLE, MPI_SUM, sim_comm);
  double RMIN = sim_RMIN(theSim);
  double RMAX = sim_RMAX(theSim);
  double ZMIN = sim_ZMIN(theSim);
  double ZMAX = sim_ZMAX(theSim);

  for (n=0;n<NUM_SCAL;++n){
    ScalarDiag_reduce[n] *= dtout/((ZMAX-ZMIN)*(RMAX*RMAX-RMIN*RMIN));
  }

  for (i=0;i<num_r_points_global;++i){
    for (n=0;n<NUM_VEC;++n){
      VectorDiag_reduce[i*NUM_VEC+n] *= dtout/(ZMAX-ZMIN);
    }
  }

  //We are doing time averaged diagnostics, so mult by delta t and add it
  //We will divide by the total delta next time we save to disk;
  for (i=0;i<num_r_points_global;++i){
    for (n=0;n<NUM_VEC;++n){
      theDiagnostics->VectorDiag[i][n] += VectorDiag_reduce[i*NUM_VEC+n] ;
    }
  }
  for (n=0;n<NUM_SCAL;++n){
    theDiagnostics->ScalarDiag[n] += ScalarDiag_reduce[n] ;
  }

  //update output time;
  theDiagnostics->toutprev = timestep_get_t(theTimeStep);

  free(ScalarDiag_temp);
  free(ScalarDiag_reduce);
  free(VectorDiag_temp);
  free(VectorDiag_reduce);
}

