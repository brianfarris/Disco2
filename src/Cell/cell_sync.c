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

//easy way to find the size of the buffer needed. Probably this could be made more efficient
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

//copy values from theCells to a buffer
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

// copy values from buffer to theCells
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

/*
void add_net_z_flux(int imin, int imax, int kmin, int kmax, struct Sim * theSim, struct Cell *** theCells){
  int NUM_Q= sim_NUM_Q(theSim);
  int count=0;
  int i,j,k,q;
  for (k=kmin;k<kmax;++k){
    for (i=imin;i<imax;++i){
      for(j = 0; j < sim_N_p(theSim,i); j++){
        theCells[k][i][j].prim[ARR] += 10.0;
        theCells[k][i][j].prim[APP] += 10.0;
        theCells[k][i][j].prim[AZZ] += 10.0;
        theCells[k][i][j].RKcons[ARR] += 10.0;
        theCells[k][i][j].RKcons[APP] += 10.0;
        theCells[k][i][j].RKcons[AZZ] += 10.0;
      }  
    }
  }

}
void subtract_net_z_flux(int imin, int imax, int kmin, int kmax, struct Sim * theSim, struct Cell *** theCells){
  int NUM_Q= sim_NUM_Q(theSim);
  int count=0;
  int i,j,k,q;
  for (k=kmin;k<kmax;++k){
    for (i=imin;i<imax;++i){
      for(j = 0; j < sim_N_p(theSim,i); j++){
        theCells[k][i][j].prim[ARR] -= 10.0;
        theCells[k][i][j].prim[APP] -= 10.0;
        theCells[k][i][j].prim[AZZ] -= 10.0;
        theCells[k][i][j].RKcons[ARR] -= 10.0;
        theCells[k][i][j].RKcons[APP] -= 10.0;
        theCells[k][i][j].RKcons[AZZ] -= 10.0;
      }  
    }
  }

}
*/


void cell_syncproc_r( struct Cell *** theCells , struct Sim *theSim,struct MPIsetup * theMPIsetup){

  //set indices for convenience
  int iNm2g = sim_N(theSim,R_DIR)-2*sim_Nghost_max(theSim,R_DIR);
  int iNmg  = sim_N(theSim,R_DIR)-sim_Nghost_max(theSim,R_DIR);
  int iN    = sim_N(theSim,R_DIR);
  int i0  = 0;
  int ig  =   sim_Nghost_min(theSim,R_DIR);
  int i2g = 2*sim_Nghost_min(theSim,R_DIR);
  int k0 = 0;
  int kN = sim_N(theSim,Z_DIR);

  int buffersize_r_hi_send,buffersize_r_hi_recv;
  double * buffer_r_hi_send;
  double * buffer_r_hi_recv;

  //if we are not at the global outer boundary, set up buffers normally
  if (!mpisetup_check_rout_bndry(theMPIsetup)){
    // find buffer sizes
    buffersize_r_hi_send = get_buffersize(iNm2g,iNmg,k0,kN,theSim);
    buffersize_r_hi_recv = get_buffersize(iNmg ,iN  ,k0,kN,theSim);
    //allocate memory for buffers
    buffer_r_hi_send = malloc(sizeof(double)*buffersize_r_hi_send);
    buffer_r_hi_recv = malloc(sizeof(double)*buffersize_r_hi_recv);
    // set buffer with data from outer radial edge of this processor
    set_buffer(iNm2g,iNmg,k0,kN,theSim,theCells,buffer_r_hi_send);

  } else{ // if we are at global outer boundary, set up a small dummy buffer
    buffersize_r_hi_send = 1;
    buffersize_r_hi_recv = 1;
    buffer_r_hi_send = malloc(sizeof(double));
    buffer_r_hi_recv = malloc(sizeof(double));
    buffer_r_hi_send[0] = 0.0;
  }

  int buffersize_r_low_send,buffersize_r_low_recv;
  double * buffer_r_low_send;
  double * buffer_r_low_recv;

  //if we are not at the global inner boundary, set up buffers normally
  if (!mpisetup_check_rin_bndry(theMPIsetup)){  
    // find buffer sizes
    buffersize_r_low_send = get_buffersize(ig,i2g,k0,kN,theSim);
    buffersize_r_low_recv = get_buffersize(i0,ig ,k0,kN,theSim);
    // allocate memory for buffers
    buffer_r_low_send = malloc(sizeof(double)*buffersize_r_low_send);
    buffer_r_low_recv = malloc(sizeof(double)*buffersize_r_low_recv);
    // set buffer with data from inner radial edge of this processor
    set_buffer(ig,i2g,k0,kN,theSim,theCells,buffer_r_low_send);
  } else{ // if we are at global inner boundary, set up a small dummy buffer
    buffersize_r_low_send = 1;
    buffersize_r_low_recv = 1;
    buffer_r_low_send = malloc(sizeof(double));
    buffer_r_low_recv = malloc(sizeof(double));
    buffer_r_low_send[0] = 0.0;
  }

  MPI_Status status;
  // send your inner buffer to the inward proc. recieve your outer buffer from the outward proc.
  MPI_Sendrecv(buffer_r_low_send,buffersize_r_low_send,MPI_DOUBLE,mpisetup_left_Proc(theMPIsetup,R_DIR),12,
      buffer_r_hi_recv,buffersize_r_hi_recv,MPI_DOUBLE,mpisetup_right_Proc(theMPIsetup,R_DIR),12,sim_comm,&status);
  // send your outer buffer to the outward proc. recieve your inner buffer from the inward proc.
  MPI_Sendrecv(buffer_r_hi_send,buffersize_r_hi_send,MPI_DOUBLE,mpisetup_right_Proc(theMPIsetup,R_DIR),13,
      buffer_r_low_recv,buffersize_r_low_recv,MPI_DOUBLE,mpisetup_left_Proc(theMPIsetup,R_DIR),13,sim_comm,&status);

  // if you are not at a global outer boundary, fill your outer ghost zones with the buffer you recieved
  if (!mpisetup_check_rout_bndry(theMPIsetup)){
    set_cells(iNmg,iN,k0,kN,theSim,theCells,buffer_r_hi_recv);
  }

  // if you are not at a global inner boundary, fill your inner ghost zones with the buffer you recieved
  if (!mpisetup_check_rin_bndry(theMPIsetup)){
    set_cells(i0,ig,k0,kN,theSim,theCells,buffer_r_low_recv);
  }

  //clean up
  free(buffer_r_low_send);
  free(buffer_r_low_recv);    
  free(buffer_r_hi_send);
  free(buffer_r_hi_recv);
}

void cell_syncproc_z( struct Cell *** theCells , struct Sim *theSim,struct MPIsetup * theMPIsetup){
  if (sim_N_global(theSim,Z_DIR)>1){
    //set indices for convenience
    int kNm2g = sim_N(theSim,Z_DIR)-2*sim_Nghost_max(theSim,Z_DIR);
    int kNmg  = sim_N(theSim,Z_DIR)-sim_Nghost_max(theSim,Z_DIR);
    int kN    = sim_N(theSim,Z_DIR);
    int k0  = 0;
    int kg  = sim_Nghost_min(theSim,Z_DIR);
    int k2g = 2*sim_Nghost_min(theSim,Z_DIR);
    int i0 = 0;
    int iN = sim_N(theSim,R_DIR);

    // find buffer sizes
    int buffersize_z_hi_send = get_buffersize(i0,iN, kNm2g,kNmg,theSim);
    int buffersize_z_hi_recv = get_buffersize(i0,iN,kNmg,kN,theSim);
    // allocate memory for buffers
    double * buffer_z_hi_send = malloc(sizeof(double)*buffersize_z_hi_send);
    double * buffer_z_hi_recv = malloc(sizeof(double)*buffersize_z_hi_recv);
    // set buffer with data from upper edge of this processor
    set_buffer(i0,iN, kNm2g,kNmg,theSim,theCells,buffer_z_hi_send);

    // find buffer sizes
    int buffersize_z_low_send = get_buffersize(i0,iN,kg, k2g,theSim);
    int buffersize_z_low_recv = get_buffersize(i0,iN,k0,kg,theSim);
    // allocate memory for buffers
    double * buffer_z_low_send = malloc(sizeof(double)*buffersize_z_low_send);
    double * buffer_z_low_recv = malloc(sizeof(double)*buffersize_z_low_recv);
    // set buffer with data from lower edge of this processor
    set_buffer(i0,iN,kg,k2g,theSim,theCells,buffer_z_low_send);

    MPI_Status status;
    // send your lower buffer to the downward proc. recieve your upper buffer from the upward proc.
    MPI_Sendrecv(buffer_z_low_send,buffersize_z_low_send,MPI_DOUBLE,mpisetup_left_Proc(theMPIsetup,Z_DIR),14,
        buffer_z_hi_recv,buffersize_z_hi_recv,MPI_DOUBLE,mpisetup_right_Proc(theMPIsetup,Z_DIR),14,sim_comm,&status);
    // send your upper buffer to the upward proc. recieve your lower buffer from the downward proc.
    MPI_Sendrecv(buffer_z_hi_send,buffersize_z_hi_send,MPI_DOUBLE,mpisetup_right_Proc(theMPIsetup,Z_DIR),15,
        buffer_z_low_recv,buffersize_z_low_recv,MPI_DOUBLE,mpisetup_left_Proc(theMPIsetup,Z_DIR),15,sim_comm,&status);

    // fill your upper ghost zones with the buffer you recieved
    set_cells(i0,iN,kNmg,kN,theSim,theCells,buffer_z_hi_recv);
    // fill your lower ghost zones with the buffer you recieved
    set_cells(i0,iN,k0,kg,theSim,theCells,buffer_z_low_recv);


    /*
       if (mpisetup_check_ztop_bndry(theMPIsetup)){
       add_net_z_flux(i0,iN,kNmg,kN,theSim,theCells);
       }

       if (mpisetup_check_zbot_bndry(theMPIsetup)){
       subtract_net_z_flux(i0,iN,kNmg,kN,theSim,theCells);
       }
       */



    //cleanup
    free(buffer_z_low_send);
    free(buffer_z_low_recv);    
    free(buffer_z_hi_send);
    free(buffer_z_hi_recv);
  }
}


