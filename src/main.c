#include <stdio.h>
#include <stdlib.h>
#include "Headers/MPIsetup.h"
#include "Headers/Cell.h"
#include "Headers/Grid.h"
#include "Headers/Face.h"
#include "Headers/GravMass.h"
#include "Headers/IO.h"
#include "Headers/TimeStep.h"
#include "Headers/header.h"
//int mpiSetup( int, char ** );

int main(int argc, char **argv) {

  // start MPI 
  // int error = mpiSetup(argc,argv);
  // if (error==1) return(0);
  struct MPIsetup * theMPIsetup = mpisetup_create(argc,argv);
  mpisetup_setprocs(theMPIsetup);
  mpisetup_cart_create(theMPIsetup);
  mpisetup_left_right(theMPIsetup);

  //grid
  struct Grid * theGrid = grid_create(theMPIsetup);
  grid_set_N_p(theGrid);
  grid_set_rz(theGrid,theMPIsetup);
  grid_set_Ncells_and_offset(theGrid,theMPIsetup);

  // gravMass
  struct GravMass *theGravMasses = gravMass_create(NP);
  gravMass_initialize(theGravMasses);
  gravMass_clean_pi(theGravMasses);

  // allocate memory for data 
  struct Cell ***theCells = cell_create(theGrid,theMPIsetup);
  cell_clean_pi(theCells,theGrid);
  // set initial data 
  cell_init(theCells,theGrid,theMPIsetup);

  //inter-processor syncs
  cell_syncproc_r(theCells,theGrid,theMPIsetup);
  cell_syncproc_z(theCells,theGrid,theMPIsetup);

  //set conserved quantities
  cell_prim2cons(theCells,theGrid);

  double dtcheck = T_MAX/NUM_CHECKPOINTS;
  double tcheck = dtcheck;

  int nfile=0;
  char filename[256];

  struct TimeStep * theTimeStep = timestep_create();
  while( timestep_get_t(theTimeStep) < T_MAX ){
    timestep_set_dt(theTimeStep,theCells,theGrid);
    cell_copy(theCells,theGrid);
    gravMass_copy(theGravMasses);
    timestep_set_RK(theTimeStep,0.0);
    timestep_substep(theTimeStep,theCells,theGrid,theGravMasses,theMPIsetup,1.0);
    timestep_set_RK(theTimeStep,0.5);
    timestep_substep(theTimeStep,theCells,theGrid,theGravMasses,theMPIsetup,0.5);
    timestep_update_Psi(theTimeStep,theCells,theGrid,theMPIsetup);
    timestep_update_t(theTimeStep); 
    if( timestep_get_t(theTimeStep)>tcheck){
      sprintf(filename,"checkpoint_%04d.h5",nfile);
      struct io *outPrims = io_create(theGrid);
      io_flattened_prim(outPrims,theCells,theGrid);
      io_hdf5_out(outPrims,theGrid,filename);
      io_destroy(outPrims,theGrid);
      tcheck += dtcheck;
      ++nfile;
    }
  }
  /*
     cell_printscreen(theCells,theGrid);
     struct io *outPrims = io_new(theGrid);
     io_flattened_prim(outPrims,theCells,theGrid);

     io_hdf5_out(outPrims,theGrid);
     io_delete(outPrims,theGrid);

     struct io *input_prims = io_new(theGrid);
     io_hdf5_in(input_prims,theGrid);
     io_unflattened_prim(outPrims,theCells,theGrid);
     io_delete(input_prims,theGrid);
     */
  //inter-processor syncs
  cell_syncproc_r(theCells,theGrid,theMPIsetup);
  cell_syncproc_z(theCells,theGrid,theMPIsetup);

  // clean up
  cell_destroy(theCells,theGrid);
  grid_destroy(theGrid);
  gravMass_destroy(theGravMasses);

  // exit MPI 
  //MPI_Barrier(grid_comm);
  //MPI_Finalize();
  mpisetup_destroy(theMPIsetup);
  return(0);
}
