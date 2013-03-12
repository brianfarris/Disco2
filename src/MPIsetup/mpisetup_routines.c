#define MPISETUP_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include "../Headers/MPIsetup.h"
#include "../Headers/header.h"

void mpisetup_setprocs(struct MPIsetup * theMPIsetup,char * filename){
  int readvar(char *,char *, int,void *);
  int err=0;
  int N_r_global,N_z_global;
  readvar( filename , "NumR" , VAR_INT  , &(N_r_global) );
  readvar( filename , "NumZ" , VAR_INT  , &(N_z_global  ));
   
  int ierr;
  int NumProcs,MyProc;
  //int dim_NumProcs[2];
  theMPIsetup->dim_NumProcs[0] = 0;
  theMPIsetup->dim_NumProcs[1] = 0;

  if (N_z_global==1){
    theMPIsetup->ndims_mpi = 1;
    theMPIsetup->dim_NumProcs[1] = 1;
    theMPIsetup->dim_MyProc[1] = 0;
  } else{
    theMPIsetup->ndims_mpi = 2;
  }
  
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &NumProcs);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &MyProc);
  //This gives me the most efficient factorization of
  //the number of processes being used into a j-by-k Mesh.
  MPI_Dims_create(NumProcs,theMPIsetup->ndims_mpi,theMPIsetup->dim_NumProcs);
  //Just in case the dimensions don't make sense; I don't know how
  //the MPI_Dims_create funcition works, so you never know.
  if( (theMPIsetup->dim_NumProcs[0])*(theMPIsetup->dim_NumProcs[1]) != NumProcs ){
    printf("Error: dimensions don't jive.\n");
    exit(1);
  }
  theMPIsetup->NumProcs=NumProcs;
  theMPIsetup->MyProc=MyProc;
}

void mpisetup_cart_create(struct MPIsetup * theMPIsetup){
  //Creating the cart!  Woohoo.
  MPI_Cart_create(MPI_COMM_WORLD,
      theMPIsetup->ndims_mpi,theMPIsetup->dim_NumProcs,theMPIsetup->wraparound,theMPIsetup->reorder,&sim_comm);
  MPI_Comm_rank(sim_comm,&(theMPIsetup->MyProc));
  MPI_Cart_coords(sim_comm,theMPIsetup->MyProc,theMPIsetup->ndims_mpi,theMPIsetup->dim_MyProc);
}

void mpisetup_left_right(struct MPIsetup * theMPIsetup){
  //Determine ranks to my left and right in all 3 dimensions.
  int next_Proc[2];
  int i;
  for( i=0 ; i<2 ; ++i ){
    next_Proc[i] = theMPIsetup->dim_MyProc[i];
  }
  for( i=0 ; i<2 ; ++i ){
    next_Proc[i] += 1;
    MPI_Cart_rank(sim_comm,next_Proc,&(theMPIsetup->right_Proc[i]));
    next_Proc[i] -= 2;
    MPI_Cart_rank(sim_comm,next_Proc,&(theMPIsetup->left_Proc[i]));
    next_Proc[i] += 1;
  }
}

