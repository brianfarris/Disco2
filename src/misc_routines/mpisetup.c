#include <stdlib.h>
#include <stdio.h>
#include "../Headers/header.h"

int mpiSetup( int argc, char **argv ){
  /* get the processor information */
  MPI_Init(&argc, &argv);
  int ierr;
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &NumProcs);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &MyProc);

  // A bunch of crap I need to create the cart. 
  int wraparound[2] = {1,1};
  dim_NumProcs[0] = 0;
  dim_NumProcs[1] = 0;
  int reorder = 1, ndims_mpi = 2;

  //This gives me the most efficient factorization of
  //the number of processes being used into a j-by-k grid.
  MPI_Dims_create(NumProcs,ndims_mpi,dim_NumProcs);

  //Just in case the dimensions don't make sense; I don't know how
  //the MPI_Dims_create funcition works, so you never know.
  if( dim_NumProcs[0]*dim_NumProcs[1] != NumProcs ){
    printf("Error: dimensions don't jive.\n");
    return(1);
  }

  //Creating the cart!  Woohoo.
  MPI_Cart_create(MPI_COMM_WORLD,ndims_mpi,dim_NumProcs,wraparound,reorder,&grid_comm);
  MPI_Comm_rank(grid_comm,&MyProc);
  MPI_Cart_coords(grid_comm,MyProc,ndims_mpi,dim_MyProc);

  //Determine ranks to my left and right in all 3 dimensions.
   int next_Proc[2];
   int i;
   for( i=0 ; i<2 ; ++i ){
      next_Proc[i] = dim_MyProc[i];
   }
   for( i=0 ; i<2 ; ++i ){
      next_Proc[i] += 1;
      MPI_Cart_rank(grid_comm,next_Proc,&(right_Proc[i]));
      next_Proc[i] -= 2;
      MPI_Cart_rank(grid_comm,next_Proc,&(left_Proc[i]));
      next_Proc[i] += 1;
   }

  return(0);
} 
