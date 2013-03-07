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

  int buffersize_r_hi_send,buffersize_r_hi_recv;
  double * buffer_r_hi_send;
  double * buffer_r_hi_recv;
  
  if (!mpisetup_check_rout_bndry(theMPIsetup)){  
    buffersize_r_hi_send = get_buffersize(grid_N_r(theGrid)-2*grid_ng(theGrid),grid_N_r(theGrid)-grid_ng(theGrid),0,grid_N_z(theGrid),theGrid);
    buffersize_r_hi_recv = get_buffersize(grid_N_r(theGrid)-grid_ng(theGrid),grid_N_r(theGrid),0,grid_N_z(theGrid),theGrid);
    buffer_r_hi_send = malloc(sizeof(double)*buffersize_r_hi_send);
    set_buffer(grid_N_r(theGrid)-2*grid_ng(theGrid),grid_N_r(theGrid)-grid_ng(theGrid),0,grid_N_z(theGrid),theGrid,theCells,buffer_r_hi_send);
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
  buffersize_r_low_send = get_buffersize(grid_ng(theGrid),2*grid_ng(theGrid),0,grid_N_z(theGrid),theGrid);
  buffersize_r_low_recv = get_buffersize(0,grid_ng(theGrid),0,grid_N_z(theGrid),theGrid);
  buffer_r_low_send = malloc(sizeof(double)*buffersize_r_low_send);
  set_buffer(grid_ng(theGrid), 2*grid_ng(theGrid),0,grid_N_z(theGrid),theGrid,theCells,buffer_r_low_send);
  buffer_r_low_recv = malloc(sizeof(double)*buffersize_r_low_recv);
  } else{
    buffersize_r_low_send = 1;
    buffersize_r_low_recv = 1;
     buffer_r_low_send = malloc(sizeof(double));
    buffer_r_low_send[0] = 0.0;
    buffer_r_low_recv = malloc(sizeof(double));
  }

  MPI_Status status;
  MPI_Sendrecv(buffer_r_low_send,buffersize_r_low_send,MPI_DOUBLE,mpisetup_left_Proc(theMPIsetup)[0],12,buffer_r_hi_recv,buffersize_r_hi_recv,MPI_DOUBLE,mpisetup_right_Proc(theMPIsetup)[0],12,grid_comm,&status);
  MPI_Sendrecv(buffer_r_hi_send,buffersize_r_hi_send,MPI_DOUBLE,mpisetup_right_Proc(theMPIsetup)[0],13,buffer_r_low_recv,buffersize_r_low_recv,MPI_DOUBLE,mpisetup_left_Proc(theMPIsetup)[0],13,grid_comm,&status);

  if (!mpisetup_check_rout_bndry(theMPIsetup)){
    set_cells(grid_N_r(theGrid)-grid_Nghost_rmax(theGrid),grid_N_r(theGrid),0,grid_N_z(theGrid),theGrid,theCells,buffer_r_hi_recv);
  }

  if (!mpisetup_check_rin_bndry(theMPIsetup)){
    set_cells(0,grid_Nghost_rmin(theGrid),0,grid_N_z(theGrid),theGrid,theCells,buffer_r_low_recv);
  }

  free(buffer_r_low_send);
  free(buffer_r_low_recv);    
  free(buffer_r_hi_send);
  free(buffer_r_hi_recv);
}

void cell_syncproc_z( struct Cell *** theCells , struct Grid *theGrid,struct MPIsetup * theMPIsetup){
  int NUM_Q = grid_NUM_Q(theGrid);

  int i,j,k,q;

  int buffersize_z_hi_send = get_buffersize(0,grid_N_r(theGrid), grid_N_z(theGrid)-2*grid_Nghost_zmax(theGrid),grid_N_z(theGrid)-grid_Nghost_zmax(theGrid),theGrid);
  int buffersize_z_hi_recv = get_buffersize(0,grid_N_r(theGrid),grid_N_z(theGrid)-grid_Nghost_zmax(theGrid),grid_N_z(theGrid),theGrid);
  double * buffer_z_hi_send = malloc(sizeof(double)*buffersize_z_hi_send);
  double * buffer_z_hi_recv = malloc(sizeof(double)*buffersize_z_hi_recv);
  set_buffer(0,grid_N_r(theGrid), grid_N_z(theGrid)-2*grid_Nghost_zmax(theGrid),grid_N_z(theGrid)-grid_Nghost_zmax(theGrid),theGrid,theCells,buffer_z_hi_send);

  int buffersize_z_low_send = get_buffersize(0,grid_N_r(theGrid),grid_Nghost_zmin(theGrid), 2*grid_Nghost_zmin(theGrid),theGrid);
  int buffersize_z_low_recv = get_buffersize(0,grid_N_r(theGrid),0,grid_Nghost_zmin(theGrid),theGrid);
  double * buffer_z_low_send = malloc(sizeof(double)*buffersize_z_low_send);
  double * buffer_z_low_recv = malloc(sizeof(double)*buffersize_z_low_recv);
  set_buffer(0,grid_N_r(theGrid), grid_Nghost_zmin(theGrid), 2*grid_Nghost_zmin(theGrid),theGrid,theCells,buffer_z_low_send);

  MPI_Status status;
  MPI_Sendrecv(buffer_z_low_send,buffersize_z_low_send,MPI_DOUBLE,mpisetup_left_Proc(theMPIsetup)[1],14,buffer_z_hi_recv,buffersize_z_hi_recv,MPI_DOUBLE,mpisetup_right_Proc(theMPIsetup)[1],14,grid_comm,&status);
  MPI_Sendrecv(buffer_z_hi_send,buffersize_z_hi_send,MPI_DOUBLE,mpisetup_right_Proc(theMPIsetup)[1],15,buffer_z_low_recv,buffersize_z_low_recv,MPI_DOUBLE,mpisetup_left_Proc(theMPIsetup)[1],15,grid_comm,&status);

  set_cells(0,grid_N_r(theGrid),grid_N_z(theGrid)-grid_Nghost_zmax(theGrid),grid_N_z(theGrid),theGrid,theCells,buffer_z_hi_recv);
  set_cells(0,grid_N_r(theGrid),0,grid_Nghost_zmin(theGrid),theGrid,theCells,buffer_z_low_recv);

  free(buffer_z_low_send);
  free(buffer_z_low_recv);    
  free(buffer_z_hi_send);
  free(buffer_z_hi_recv);
}


