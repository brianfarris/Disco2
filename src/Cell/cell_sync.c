#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Sim.h"
#include "../Headers/Face.h"
#include "../Headers/GravMass.h"
#include "../Headers/MPIsetup.h"
#include "../Headers/header.h"

int get_buffersize(int imin, int imax, int kmin, int kmax, struct Sim * theSim){
  int NUM_Q= sim_NUM_Q(theSim);
  int size=0;
  int i,j,k;
  for (k=kmin;k<kmax;++k){
    for (i=imin;i<imax;++i){
      for(j = 0; j < sim_N_p(theSim,i); j++){
        size++;
      }
    }
  }
  size *=(2*NUM_Q+1);
  return(size);
}

void set_buffer(int imin, int imax, int kmin, int kmax, struct Sim * theSim, struct Cell *** theCells, double * buffer){
  int NUM_Q= sim_NUM_Q(theSim);
  int count=0;
  int i,j,k,q;
  for (k=kmin;k<kmax;++k){
    for (i=imin;i<imax;++i){
      for(j = 0; j < sim_N_p(theSim,i); j++){
        for (q=0;q<NUM_Q;++q){
          buffer[count] = theCells[k][i][j].prim[q];
          ++count;
        }
        for (q=0;q<NUM_Q;++q){
          buffer[count] = theCells[k][i][j].RKcons[q];
          ++count;
        }
        buffer[count] = theCells[k][i][j].tiph;
        ++count;
      }
    }
  }
}

void set_cells(int imin, int imax, int kmin, int kmax, struct Sim * theSim, struct Cell *** theCells, double * buffer){
  int NUM_Q= sim_NUM_Q(theSim);
  int count=0;
  int i,j,k,q;
  for (k=kmin;k<kmax;++k){
    for (i=imin;i<imax;++i){
      for(j = 0; j < sim_N_p(theSim,i); j++){
        for (q=0;q<NUM_Q;++q){
          theCells[k][i][j].prim[q] = buffer[count];
          ++count;
        }  
        for (q=0;q<NUM_Q;++q){
          theCells[k][i][j].RKcons[q] = buffer[count];
          ++count;
        }
        theCells[k][i][j].tiph = buffer[count];
        ++count;
      }
    }
  }
}

void cell_syncproc_r( struct Cell *** theCells , struct Sim *theSim,struct MPIsetup * theMPIsetup){

  int buffersize_r_hi_send,buffersize_r_hi_recv;
  double * buffer_r_hi_send;
  double * buffer_r_hi_recv;
  
  if (!mpisetup_check_rout_bndry(theMPIsetup)){  
    buffersize_r_hi_send = get_buffersize(sim_N(theSim,R_DIR)-2*sim_ng(theSim),sim_N(theSim,R_DIR)-sim_ng(theSim),0,sim_N(theSim,Z_DIR),theSim);
    buffersize_r_hi_recv = get_buffersize(sim_N(theSim,R_DIR)-sim_ng(theSim),sim_N(theSim,R_DIR),0,sim_N(theSim,Z_DIR),theSim);
    buffer_r_hi_send = malloc(sizeof(double)*buffersize_r_hi_send);
    set_buffer(sim_N(theSim,R_DIR)-2*sim_ng(theSim),sim_N(theSim,R_DIR)-sim_ng(theSim),0,sim_N(theSim,Z_DIR),theSim,theCells,buffer_r_hi_send);
    buffer_r_hi_recv = malloc(sizeof(double)*buffersize_r_hi_recv);
  } else{
    buffersize_r_hi_send = 1;
    buffersize_r_hi_recv = 1;
    buffer_r_hi_send = malloc(sizeof(double));
    buffer_r_hi_send[0] = 0.0;
    buffer_r_hi_recv = malloc(sizeof(double));
  }

  int buffersize_r_low_send,buffersize_r_low_recv;
  double * buffer_r_low_send;
  double * buffer_r_low_recv;

  if (!mpisetup_check_rin_bndry(theMPIsetup)){  
  buffersize_r_low_send = get_buffersize(sim_ng(theSim),2*sim_ng(theSim),0,sim_N(theSim,Z_DIR),theSim);
  buffersize_r_low_recv = get_buffersize(0,sim_ng(theSim),0,sim_N(theSim,Z_DIR),theSim);
  buffer_r_low_send = malloc(sizeof(double)*buffersize_r_low_send);
  set_buffer(sim_ng(theSim), 2*sim_ng(theSim),0,sim_N(theSim,Z_DIR),theSim,theCells,buffer_r_low_send);
  buffer_r_low_recv = malloc(sizeof(double)*buffersize_r_low_recv);
  } else{
    buffersize_r_low_send = 1;
    buffersize_r_low_recv = 1;
     buffer_r_low_send = malloc(sizeof(double));
    buffer_r_low_send[0] = 0.0;
    buffer_r_low_recv = malloc(sizeof(double));
  }

  MPI_Status status;
  MPI_Sendrecv(buffer_r_low_send,buffersize_r_low_send,MPI_DOUBLE,mpisetup_left_Proc(theMPIsetup)[0],12,buffer_r_hi_recv,buffersize_r_hi_recv,MPI_DOUBLE,mpisetup_right_Proc(theMPIsetup)[0],12,sim_comm,&status);
  MPI_Sendrecv(buffer_r_hi_send,buffersize_r_hi_send,MPI_DOUBLE,mpisetup_right_Proc(theMPIsetup)[0],13,buffer_r_low_recv,buffersize_r_low_recv,MPI_DOUBLE,mpisetup_left_Proc(theMPIsetup)[0],13,sim_comm,&status);

  if (!mpisetup_check_rout_bndry(theMPIsetup)){
    set_cells(sim_N(theSim,R_DIR)-sim_Nghost_max(theSim,R_DIR),sim_N(theSim,R_DIR),0,sim_N(theSim,Z_DIR),theSim,theCells,buffer_r_hi_recv);
  }

  if (!mpisetup_check_rin_bndry(theMPIsetup)){
    set_cells(0,sim_Nghost_min(theSim,R_DIR),0,sim_N(theSim,Z_DIR),theSim,theCells,buffer_r_low_recv);
  }

  free(buffer_r_low_send);
  free(buffer_r_low_recv);    
  free(buffer_r_hi_send);
  free(buffer_r_hi_recv);
}

void cell_syncproc_z( struct Cell *** theCells , struct Sim *theSim,struct MPIsetup * theMPIsetup){
  if (sim_N_global(theSim,Z_DIR)>1){
    int NUM_Q = sim_NUM_Q(theSim);

    int i,j,k,q;

    int buffersize_z_hi_send = get_buffersize(0,sim_N(theSim,R_DIR), sim_N(theSim,Z_DIR)-2*sim_Nghost_max(theSim,Z_DIR),sim_N(theSim,Z_DIR)-sim_Nghost_max(theSim,Z_DIR),theSim);
    int buffersize_z_hi_recv = get_buffersize(0,sim_N(theSim,R_DIR),sim_N(theSim,Z_DIR)-sim_Nghost_max(theSim,Z_DIR),sim_N(theSim,Z_DIR),theSim);
    double * buffer_z_hi_send = malloc(sizeof(double)*buffersize_z_hi_send);
    double * buffer_z_hi_recv = malloc(sizeof(double)*buffersize_z_hi_recv);
    set_buffer(0,sim_N(theSim,R_DIR), sim_N(theSim,Z_DIR)-2*sim_Nghost_max(theSim,Z_DIR),sim_N(theSim,Z_DIR)-sim_Nghost_max(theSim,Z_DIR),theSim,theCells,buffer_z_hi_send);

    int buffersize_z_low_send = get_buffersize(0,sim_N(theSim,R_DIR),sim_Nghost_min(theSim,Z_DIR), 2*sim_Nghost_min(theSim,Z_DIR),theSim);
    int buffersize_z_low_recv = get_buffersize(0,sim_N(theSim,R_DIR),0,sim_Nghost_min(theSim,Z_DIR),theSim);
    double * buffer_z_low_send = malloc(sizeof(double)*buffersize_z_low_send);
    double * buffer_z_low_recv = malloc(sizeof(double)*buffersize_z_low_recv);
    set_buffer(0,sim_N(theSim,R_DIR), sim_Nghost_min(theSim,Z_DIR), 2*sim_Nghost_min(theSim,Z_DIR),theSim,theCells,buffer_z_low_send);

    MPI_Status status;
    MPI_Sendrecv(buffer_z_low_send,buffersize_z_low_send,MPI_DOUBLE,mpisetup_left_Proc(theMPIsetup)[1],14,buffer_z_hi_recv,buffersize_z_hi_recv,MPI_DOUBLE,mpisetup_right_Proc(theMPIsetup)[1],14,sim_comm,&status);
    MPI_Sendrecv(buffer_z_hi_send,buffersize_z_hi_send,MPI_DOUBLE,mpisetup_right_Proc(theMPIsetup)[1],15,buffer_z_low_recv,buffersize_z_low_recv,MPI_DOUBLE,mpisetup_left_Proc(theMPIsetup)[1],15,sim_comm,&status);

    set_cells(0,sim_N(theSim,R_DIR),sim_N(theSim,Z_DIR)-sim_Nghost_max(theSim,Z_DIR),sim_N(theSim,Z_DIR),theSim,theCells,buffer_z_hi_recv);
    set_cells(0,sim_N(theSim,R_DIR),0,sim_Nghost_min(theSim,Z_DIR),theSim,theCells,buffer_z_low_recv);

    free(buffer_z_low_send);
    free(buffer_z_low_recv);    
    free(buffer_z_hi_send);
    free(buffer_z_hi_recv);
  }
}


