#define DIAGNOSTICS_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Sim.h"
#include "../Headers/Diagnostics.h"
#include "../Headers/MPIsetup.h"
#include "../Headers/TimeStep.h"
#include "../Headers/header.h"

struct Diagnostics *diagnostics_create(struct Sim * theSim, struct TimeStep * theTimeStep, struct MPIsetup * theMPIsetup) {
  struct Diagnostics * theDiagnostics = (struct Diagnostics *) malloc(sizeof(struct Diagnostics));

  theDiagnostics->NUM_DIAG = 16; //DD change this if add new output in diagnostics_set (was 16)
  //int NUM_TST  = 1;

  int * N_p_global_temp = malloc(sizeof(int)*sim_N_global(theSim,R_DIR));
  theDiagnostics->N_p_global = malloc(sizeof(int)*sim_N_global(theSim,R_DIR));
  int i,j,k;
  for (i=0;i<sim_N_global(theSim,R_DIR);i++){
    N_p_global_temp[i] = 0;
  }

  int imin = sim_Nghost_min(theSim,R_DIR);
  int imax = sim_N(theSim,R_DIR)-sim_Nghost_max(theSim,R_DIR);
  int kmin = sim_Nghost_min(theSim,Z_DIR);
  int kmax = sim_N(theSim,Z_DIR)-sim_Nghost_max(theSim,Z_DIR);


  // Stuff that involves counting cells  
  int Ncells_eq_global;
  int Ncells_eq=0;
  for (k=kmin;k<kmax;++k){
    double zp = sim_FacePos(theSim,k,Z_DIR);
    double zm = sim_FacePos(theSim,k-1,Z_DIR);
    double z = 0.5*(zm+zp);
    double dz = zp-zm;
    for (i=imin;i<imax;++i){
      double rp = sim_FacePos(theSim,i,R_DIR);
      double rm = sim_FacePos(theSim,i-1,R_DIR);
      double r = 0.5*(rm+rp);
      for (j=0;j<sim_N_p(theSim,i);++j){
        if ((fabs(zp)<0.0000001)||(fabs(z)<0.0000001)){
          ++Ncells_eq; 
        }
      }
    }
  }

  int *Ncells_eq_arr = malloc(sizeof(int) * mpisetup_NumProcs(theMPIsetup));
  MPI_Allgather(&Ncells_eq, 1, MPI_INT, Ncells_eq_arr, 1, MPI_INT, MPI_COMM_WORLD);
  int offset_eq=0;
  for (i=0;i<mpisetup_MyProc(theMPIsetup);i++){
    offset_eq += Ncells_eq_arr[i];
  }
  free(Ncells_eq_arr);

  MPI_Allreduce( &Ncells_eq , &Ncells_eq_global , 1 , MPI_INT , MPI_SUM , sim_comm );

  theDiagnostics->offset_eq = offset_eq;
  theDiagnostics->N_eq_cells = Ncells_eq_global;

  for (i=imin;i<imax;i++){
    N_p_global_temp[sim_N0(theSim,R_DIR)+i-imin] = sim_N_p(theSim,i)*mpisetup_check_ztop_bndry(theMPIsetup);
  }
  MPI_Allreduce( N_p_global_temp,theDiagnostics->N_p_global , sim_N_global(theSim,R_DIR), MPI_INT, MPI_SUM, sim_comm);
  free(N_p_global_temp);

  theDiagnostics->EquatDiag = malloc(sizeof(double *) * theDiagnostics->N_eq_cells);
  theDiagnostics->EquatDiag[0] = malloc(sizeof(double) * theDiagnostics->N_eq_cells*(theDiagnostics->NUM_DIAG+2));

  int position = 0;
  for (i=0;i<sim_N_global(theSim,R_DIR);i++){
    for (j=0;j<theDiagnostics->N_p_global[i];j++){
      theDiagnostics->EquatDiag[position] = theDiagnostics->EquatDiag[0] + position*(theDiagnostics->NUM_DIAG+2);
      ++position;
    }
  }

  theDiagnostics->VectorDiag = malloc(sizeof(double *) * sim_N_global(theSim,R_DIR));
  theDiagnostics->VectorDiag[0] = malloc(sizeof(double) * sim_N_global(theSim,R_DIR) * (theDiagnostics->NUM_DIAG+1));
  int n;
  for (i=0;i<sim_N_global(theSim,R_DIR);++i){
    theDiagnostics->VectorDiag[i] = theDiagnostics->VectorDiag[0]+i*(theDiagnostics->NUM_DIAG+1);
  }

  theDiagnostics->ScalarDiag = malloc(sizeof(double)*theDiagnostics->NUM_DIAG);

  //DD FOR TRQ TST -----
  //theDiagnostics->TrVec = malloc(sizeof(double *) * sim_N_global(theSim,R_DIR));
  //theDiagnostics->TrVec[0] = malloc(sizeof(double) * sim_N_global(theSim,R_DIR) * (NUM_TST+1));
  //for (i=0;i<sim_N_global(theSim,R_DIR);++i){
  //  theDiagnostics->TrVec[i] = theDiagnostics->TrVec[0]+i*(NUM_TST+1);
  //}

  //theDiagnostics->TrScal = malloc(sizeof(double)*NUM_TST);
  //------------------------


  position = 0;
  for (i=0;i<sim_N_global(theSim,R_DIR);++i){
    for (j=0;j<theDiagnostics->N_p_global[i];++j){
      for (n=0;n<(2+theDiagnostics->NUM_DIAG);++n){
        theDiagnostics->EquatDiag[position][n]=0.0;
      }
      ++position;
    }
  }
  for (i=0;i<sim_N_global(theSim,R_DIR);++i){
    for (n=0;n<(1+theDiagnostics->NUM_DIAG);++n){
      theDiagnostics->VectorDiag[i][n]=0.0;
    }
  }
  for (n=0;n<theDiagnostics->NUM_DIAG;++n){
    theDiagnostics->ScalarDiag[n]=0.0;
  }

  theDiagnostics->toutprev = timestep_get_t(theTimeStep); 
  theDiagnostics->toutprev_dump = timestep_get_t(theTimeStep);

  // for determining when to measure diagnostics
  theDiagnostics->dtdiag_measure = sim_get_T_MAX(theSim)/sim_NUM_DIAG_MEASURE(theSim);
  theDiagnostics->tdiag_measure = theDiagnostics->dtdiag_measure;
  // for determining when to dump diagnostics
  theDiagnostics->dtdiag_dump = sim_get_T_MAX(theSim)/sim_NUM_DIAG_DUMP(theSim);
  theDiagnostics->tdiag_dump = theDiagnostics->dtdiag_dump;

  while( theDiagnostics->tdiag_measure< timestep_get_t(theTimeStep)){
    theDiagnostics->tdiag_measure+=theDiagnostics->dtdiag_measure;
  }
  while( theDiagnostics->tdiag_dump< timestep_get_t(theTimeStep)){
    theDiagnostics->tdiag_dump+=theDiagnostics->dtdiag_dump;
  }

  return theDiagnostics;
}

void diagnostics_destroy(struct Diagnostics *theDiagnostics,struct Sim * theSim) {
  free(theDiagnostics->ScalarDiag);
  free(theDiagnostics->VectorDiag[0]);
  free(theDiagnostics->VectorDiag);
  free(theDiagnostics->EquatDiag[0]);
  free(theDiagnostics->EquatDiag);
  free(theDiagnostics->N_p_global);
  free(theDiagnostics);
  //DD
  //free(theDiagnostics->TrScal);
  //free(theDiagnostics->TrVec[0]);
  //free(theDiagnostics->TrVec);
}


