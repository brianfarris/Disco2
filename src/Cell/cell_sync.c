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

void cell_syncproc_r( struct Cell *** theCells , struct Grid *theGrid,struct MPIsetup * theMPIsetup){
  int N_r_withghost = grid_N_r(theGrid)+grid_Nghost_rmin(theGrid)+grid_Nghost_rmax(theGrid);
  int N_z_withghost = grid_N_z(theGrid)+grid_Nghost_zmin(theGrid)+grid_Nghost_zmax(theGrid);

  int i,j,k,q;

  int buffersize_r_out_send=0;
  int buffersize_r_out_recv=0;
  double *buffer_r_out_send;
  double *buffer_r_out_recv;

  int count=0;
  for (k=0;k<N_z_withghost;++k){
    for (i= N_r_withghost-2*ng;i< N_r_withghost-ng;++i){
      for(j = 0; j < grid_N_p(theGrid,i); j++){
        buffersize_r_out_send++;
      }
    }
  }
  buffersize_r_out_send *= (2*NUM_Q+1);

  for (k=0;k<N_z_withghost;++k){
    for (i= N_r_withghost-ng;i< N_r_withghost;++i){
      for(j = 0; j < grid_N_p(theGrid,i); j++){
        buffersize_r_out_recv++;
      }
    }
  }
  buffersize_r_out_recv *= (2*NUM_Q+1);

  buffer_r_out_send = malloc(sizeof(double)*buffersize_r_out_send);
  buffer_r_out_recv = malloc(sizeof(double)*buffersize_r_out_recv);
  
  count=0;
  for (k=0;k<N_z_withghost;++k){
    for (i= N_r_withghost-2*ng;i< N_r_withghost-ng;++i){
      for(j = 0; j < grid_N_p(theGrid,i); j++){
        for (q=0;q<NUM_Q;++q){
          buffer_r_out_send[count] = theCells[k][i][j].prim[q];
          ++count;
        }
        for (q=0;q<NUM_Q;++q){
          buffer_r_out_send[count] = theCells[k][i][j].RKcons[q];
          ++count;
        }
        buffer_r_out_send[count] = theCells[k][i][j].tiph;
        ++count;
      }
    }
  }

  int buffersize_r_in_send=0;
  int buffersize_r_in_recv=0;
  double *buffer_r_in_send;
  double *buffer_r_in_recv;

  for (k=0;k<N_z_withghost;++k){
    for (i= ng;i< 2*ng;++i){
      for(j = 0; j < grid_N_p(theGrid,i); j++){
        buffersize_r_in_send++;
      }
    }
  }
  buffersize_r_in_send *= (2*NUM_Q+1);
  for (k=0;k<N_z_withghost;++k){
    for (i=0;i<ng;++i){
      for(j = 0; j < grid_N_p(theGrid,i); j++){
        buffersize_r_in_recv++;
      }
    }
  }
  buffersize_r_in_recv *= (2*NUM_Q+1);

  buffer_r_in_send = malloc(sizeof(double)*buffersize_r_in_send);
  buffer_r_in_recv = malloc(sizeof(double)*buffersize_r_in_recv);
  count=0;
  for (k=0;k<N_z_withghost;++k){
    for (i= ng;i< 2*ng;++i){
      for(j = 0; j < grid_N_p(theGrid,i); j++){
        for (q=0;q<NUM_Q;++q){
          buffer_r_in_send[count] = theCells[k][i][j].prim[q];
          ++count;
        }
        for (q=0;q<NUM_Q;++q){
          buffer_r_in_send[count] = theCells[k][i][j].RKcons[q];
          ++count;
        }
        buffer_r_in_send[count] = theCells[k][i][j].tiph;
        ++count;
      }
    }
  }
  MPI_Status status;

  MPI_Sendrecv(buffer_r_in_send,buffersize_r_in_send,MPI_DOUBLE,mpisetup_left_Proc(theMPIsetup)[0],12,buffer_r_out_recv,buffersize_r_out_recv,MPI_DOUBLE,mpisetup_right_Proc(theMPIsetup)[0],12,grid_comm,&status);
  MPI_Sendrecv(buffer_r_out_send,buffersize_r_out_send,MPI_DOUBLE,mpisetup_right_Proc(theMPIsetup)[0],13,buffer_r_in_recv,buffersize_r_in_recv,MPI_DOUBLE,mpisetup_left_Proc(theMPIsetup)[0],13,grid_comm,&status);


  //if (dim_MyProc[0]!=(dim_NumProcs[0]-1)){
  if (!mpisetup_check_rout_bndry(theMPIsetup)){
    count=0;
    for (k=0;k<N_z_withghost;++k){
      for (i= N_r_withghost-grid_Nghost_rmax(theGrid);i< N_r_withghost;++i){
        for(j = 0; j < grid_N_p(theGrid,i); j++){
          for (q=0;q<NUM_Q;++q){
            theCells[k][i][j].prim[q] = buffer_r_out_recv[count];
            ++count;
          }
          for (q=0;q<NUM_Q;++q){
            theCells[k][i][j].RKcons[q] = buffer_r_out_recv[count];
            ++count;
          }
          theCells[k][i][j].tiph = buffer_r_out_recv[count];
          ++count;
        }
      }
    }
  }

  if (!mpisetup_check_rin_bndry(theMPIsetup)){
    count=0;
    for (k=0;k<N_z_withghost;++k){
      for (i=0;i<grid_Nghost_rmin(theGrid);++i){
        for(j = 0; j < grid_N_p(theGrid,i); j++){
          for (q=0;q<NUM_Q;++q){
            theCells[k][i][j].prim[q] = buffer_r_in_recv[count];
            ++count;
          }
          for (q=0;q<NUM_Q;++q){
            theCells[k][i][j].RKcons[q] = buffer_r_in_recv[count];
            ++count;
          }
          theCells[k][i][j].tiph = buffer_r_out_recv[count];
          ++count;
        }
      }
    }
  }
  free(buffer_r_in_send);
  free(buffer_r_in_recv);    
  free(buffer_r_out_send);
  free(buffer_r_out_recv);
}


void cell_syncproc_z( struct Cell *** theCells , struct Grid *theGrid,struct MPIsetup * theMPIsetup){
  int N_r_withghost = grid_N_r(theGrid)+grid_Nghost_rmin(theGrid)+grid_Nghost_rmax(theGrid);
  int N_z_withghost = grid_N_z(theGrid)+grid_Nghost_zmin(theGrid)+grid_Nghost_zmax(theGrid);

  int i,j,k,q;

  int buffersize_z_top_send=0;
  int buffersize_z_top_recv=0;
  double *buffer_z_top_send;
  double *buffer_z_top_recv;

  int count=0;
  for (k= N_z_withghost-2*grid_Nghost_zmax(theGrid);k< N_z_withghost-grid_Nghost_zmax(theGrid);++k){
    for (i=0;i<N_r_withghost;++i){
      for(j = 0; j < grid_N_p(theGrid,i); j++){
        buffersize_z_top_send++;
      }
    }
  }
  buffersize_z_top_send *= (2*NUM_Q+1);

  for (k= N_z_withghost-grid_Nghost_zmax(theGrid);k< N_z_withghost;++k){
    for (i=0;i<N_r_withghost;++i){
      for(j = 0; j < grid_N_p(theGrid,i); j++){
        buffersize_z_top_recv++;
      }
    }
  }
  buffersize_z_top_recv *= (2*NUM_Q+1);

  buffer_z_top_send = malloc(sizeof(double)*buffersize_z_top_send);
  buffer_z_top_recv = malloc(sizeof(double)*buffersize_z_top_recv);
  count=0;
  for (k= N_z_withghost-2*grid_Nghost_zmax(theGrid);k< N_z_withghost-grid_Nghost_zmax(theGrid);++k){
    for (i=0;i<N_r_withghost;++i){
      for(j = 0; j < grid_N_p(theGrid,i); j++){
        for (q=0;q<NUM_Q;++q){
          buffer_z_top_send[count] = theCells[k][i][j].prim[q];
          ++count;
        }
        for (q=0;q<NUM_Q;++q){
          buffer_z_top_send[count] = theCells[k][i][j].RKcons[q];
          ++count;
        }
        buffer_z_top_send[count] = theCells[k][i][j].tiph;
        ++count;
      }
    }
  }

  int buffersize_z_bot_send=0;
  int buffersize_z_bot_recv=0;
  double *buffer_z_bot_send;
  double *buffer_z_bot_recv;

  for (k= grid_Nghost_zmin(theGrid);k< 2*grid_Nghost_zmin(theGrid);++k){
    for (i=0;i<N_r_withghost;++i){
      for(j = 0; j < grid_N_p(theGrid,i); j++){
        buffersize_z_bot_send++;
      }
    }
  }
  buffersize_z_bot_send *= (2*NUM_Q+1);
  for (k=0;k<grid_Nghost_zmin(theGrid);++k){
    for (i=0;i<N_r_withghost;++i){
      for(j = 0; j < grid_N_p(theGrid,i); j++){
        buffersize_z_bot_recv++;
      }
    }
  }
  buffersize_z_bot_recv *= (2*NUM_Q+1);

  buffer_z_bot_send = malloc(sizeof(double)*buffersize_z_bot_send);
  buffer_z_bot_recv = malloc(sizeof(double)*buffersize_z_bot_recv);
  count=0;
  for (k= grid_Nghost_zmin(theGrid);k< 2*grid_Nghost_zmin(theGrid);++k){
    for (i=0;i<N_r_withghost;++i){
      for(j = 0; j < grid_N_p(theGrid,i); j++){
        for (q=0;q<NUM_Q;++q){
          buffer_z_bot_send[count] = theCells[k][i][j].prim[q];
          ++count;
        }
        for (q=0;q<NUM_Q;++q){
          buffer_z_bot_send[count] = theCells[k][i][j].RKcons[q];
          ++count;
        }
        buffer_z_bot_send[count] = theCells[k][i][j].tiph;
        ++count;
      }
    }
  }
  MPI_Status status;

  MPI_Sendrecv(buffer_z_bot_send,buffersize_z_bot_send,MPI_DOUBLE,mpisetup_left_Proc(theMPIsetup)[1],14,buffer_z_top_recv,buffersize_z_top_recv,MPI_DOUBLE,mpisetup_right_Proc(theMPIsetup)[1],14,grid_comm,&status);
  MPI_Sendrecv(buffer_z_top_send,buffersize_z_top_send,MPI_DOUBLE,mpisetup_right_Proc(theMPIsetup)[1],15,buffer_z_bot_recv,buffersize_z_bot_recv,MPI_DOUBLE,mpisetup_left_Proc(theMPIsetup)[1],15,grid_comm,&status);

  count=0;
  for (k= N_z_withghost-grid_Nghost_zmax(theGrid);k< N_z_withghost;++k){
    for (i=0;i<N_r_withghost;++i){
      for(j = 0; j < grid_N_p(theGrid,i); j++){
        for (q=0;q<NUM_Q;++q){
          theCells[k][i][j].prim[q] = buffer_z_top_recv[count];
          ++count;
        }
        for (q=0;q<NUM_Q;++q){
          theCells[k][i][j].RKcons[q] = buffer_z_top_recv[count];
          ++count;
        }
        theCells[k][i][j].tiph = buffer_z_top_recv[count];
        ++count;
      }
    }
  }

  count=0;
  for (k=0;k<grid_Nghost_zmin(theGrid);++k){
    for (i=0;i<N_r_withghost;++i){
      for(j = 0; j < grid_N_p(theGrid,i); j++){
        for (q=0;q<NUM_Q;++q){
          theCells[k][i][j].prim[q] = buffer_z_bot_recv[count];
          ++count;
        }
        for (q=0;q<NUM_Q;++q){
          theCells[k][i][j].RKcons[q] = buffer_z_bot_recv[count];
          ++count;
        }
        theCells[k][i][j].tiph = buffer_z_top_recv[count];
        ++count;
      }
    }
  }
  free(buffer_z_bot_send);
  free(buffer_z_bot_recv);    
  free(buffer_z_top_send);
  free(buffer_z_top_recv);
}


