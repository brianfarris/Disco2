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

int get_buffersize(int imin, int imax, int kmin, int kmax, struct Grid * theGrid){
  int NUM_Q= grid_NUM_Q(theGrid);
  int size=0;
  int i,j,k;
  for (k=kmin;k<kmax;++k){
    for (i=imin;i<imax;++i){
      for(j = 0; j < grid_N_p(theGrid,i); j++){
        size++;
      }
    }
  }
  size *=(2*NUM_Q+1);
  return(size);
}

void set_buffer(int imin, int imax, int kmin, int kmax, struct Grid * theGrid, struct Cell *** theCells, double * buffer){
  int NUM_Q= grid_NUM_Q(theGrid);
  int count=0;
  int i,j,k,q;
  for (k=kmin;k<kmax;++k){
    for (i=imin;i<imax;++i){
      for(j = 0; j < grid_N_p(theGrid,i); j++){
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

void set_cells(int imin, int imax, int kmin, int kmax, struct Grid * theGrid, struct Cell *** theCells, double * buffer){
  int NUM_Q= grid_NUM_Q(theGrid);
  int count=0;
  int i,j,k,q;
  for (k=kmin;k<kmax;++k){
    for (i=imin;i<imax;++i){
      for(j = 0; j < grid_N_p(theGrid,i); j++){
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

void cell_syncproc_r( struct Cell *** theCells , struct Grid *theGrid,struct MPIsetup * theMPIsetup){
  int NUM_Q = grid_NUM_Q(theGrid);

  int i,j,k,q;

  int buffersize_r_out_send,buffersize_r_out_recv;
  double * buffer_r_out_send;
  double * buffer_r_out_recv;
  
  if (!mpisetup_check_rout_bndry(theMPIsetup)){  
    buffersize_r_out_send = get_buffersize(grid_N_r(theGrid)-2*grid_ng(theGrid),grid_N_r(theGrid)-grid_ng(theGrid),0,grid_N_z(theGrid),theGrid);
    buffersize_r_out_recv = get_buffersize(grid_N_r(theGrid)-grid_ng(theGrid),grid_N_r(theGrid),0,grid_N_z(theGrid),theGrid);
    buffer_r_out_send = malloc(sizeof(double)*buffersize_r_out_send);
    set_buffer(grid_N_r(theGrid)-2*grid_ng(theGrid),grid_N_r(theGrid)-grid_ng(theGrid),0,grid_N_z(theGrid),theGrid,theCells,buffer_r_out_send);
    buffer_r_out_recv = malloc(sizeof(double)*buffersize_r_out_recv);
  } else{
    buffersize_r_out_send = 1;
    buffersize_r_out_recv = 1;
    buffer_r_out_send = malloc(sizeof(double));
    buffer_r_out_send[0] = 0.0;
    buffer_r_out_recv = malloc(sizeof(double));
  }

  int buffersize_r_in_send,buffersize_r_in_recv;
  double * buffer_r_in_send;
  double * buffer_r_in_recv;

  if (!mpisetup_check_rin_bndry(theMPIsetup)){  
  buffersize_r_in_send = get_buffersize(grid_ng(theGrid),2*grid_ng(theGrid),0,grid_N_z(theGrid),theGrid);
  buffersize_r_in_recv = get_buffersize(0,grid_ng(theGrid),0,grid_N_z(theGrid),theGrid);
  buffer_r_in_send = malloc(sizeof(double)*buffersize_r_in_send);
  set_buffer(grid_ng(theGrid), 2*grid_ng(theGrid),0,grid_N_z(theGrid),theGrid,theCells,buffer_r_in_send);
  buffer_r_in_recv = malloc(sizeof(double)*buffersize_r_in_recv);
  } else{
    buffersize_r_in_send = 1;
    buffersize_r_in_recv = 1;
     buffer_r_in_send = malloc(sizeof(double));
    buffer_r_in_send[0] = 0.0;
    buffer_r_in_recv = malloc(sizeof(double));
  }



  MPI_Status status;

  MPI_Sendrecv(buffer_r_in_send,buffersize_r_in_send,MPI_DOUBLE,mpisetup_left_Proc(theMPIsetup)[0],12,buffer_r_out_recv,buffersize_r_out_recv,MPI_DOUBLE,mpisetup_right_Proc(theMPIsetup)[0],12,grid_comm,&status);
  MPI_Sendrecv(buffer_r_out_send,buffersize_r_out_send,MPI_DOUBLE,mpisetup_right_Proc(theMPIsetup)[0],13,buffer_r_in_recv,buffersize_r_in_recv,MPI_DOUBLE,mpisetup_left_Proc(theMPIsetup)[0],13,grid_comm,&status);


  if (!mpisetup_check_rout_bndry(theMPIsetup)){
    set_cells(grid_N_r(theGrid)-grid_Nghost_rmax(theGrid),grid_N_r(theGrid),0,grid_N_z(theGrid),theGrid,theCells,buffer_r_out_recv);
  }

  if (!mpisetup_check_rin_bndry(theMPIsetup)){
    set_cells(0,grid_Nghost_rmin(theGrid),0,grid_N_z(theGrid),theGrid,theCells,buffer_r_in_recv);
  }

  free(buffer_r_in_send);
  free(buffer_r_in_recv);    
  free(buffer_r_out_send);
  free(buffer_r_out_recv);
}


void cell_syncproc_z( struct Cell *** theCells , struct Grid *theGrid,struct MPIsetup * theMPIsetup){
  int NUM_Q = grid_NUM_Q(theGrid);

  int i,j,k,q;

  int buffersize_z_top_send=0;
  int buffersize_z_top_recv=0;
  double *buffer_z_top_send;
  double *buffer_z_top_recv;

  int count=0;
  for (k= grid_N_z(theGrid)-2*grid_Nghost_zmax(theGrid);k< grid_N_z(theGrid)-grid_Nghost_zmax(theGrid);++k){
    for (i=0;i<grid_N_r(theGrid);++i){
      for(j = 0; j < grid_N_p(theGrid,i); j++){
        buffersize_z_top_send++;
      }
    }
  }
  buffersize_z_top_send *= (2*NUM_Q+1);

  for (k= grid_N_z(theGrid)-grid_Nghost_zmax(theGrid);k< grid_N_z(theGrid);++k){
    for (i=0;i<grid_N_r(theGrid);++i){
      for(j = 0; j < grid_N_p(theGrid,i); j++){
        buffersize_z_top_recv++;
      }
    }
  }
  buffersize_z_top_recv *= (2*NUM_Q+1);

  buffer_z_top_send = malloc(sizeof(double)*buffersize_z_top_send);
  buffer_z_top_recv = malloc(sizeof(double)*buffersize_z_top_recv);
  count=0;
  for (k= grid_N_z(theGrid)-2*grid_Nghost_zmax(theGrid);k< grid_N_z(theGrid)-grid_Nghost_zmax(theGrid);++k){
    for (i=0;i<grid_N_r(theGrid);++i){
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
    for (i=0;i<grid_N_r(theGrid);++i){
      for(j = 0; j < grid_N_p(theGrid,i); j++){
        buffersize_z_bot_send++;
      }
    }
  }
  buffersize_z_bot_send *= (2*NUM_Q+1);
  for (k=0;k<grid_Nghost_zmin(theGrid);++k){
    for (i=0;i<grid_N_r(theGrid);++i){
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
    for (i=0;i<grid_N_r(theGrid);++i){
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
  for (k= grid_N_z(theGrid)-grid_Nghost_zmax(theGrid);k< grid_N_z(theGrid);++k){
    for (i=0;i<grid_N_r(theGrid);++i){
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
    for (i=0;i<grid_N_r(theGrid);++i){
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


