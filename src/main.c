#include <stdio.h>
#include <stdlib.h>
#include "Headers/mpisetup.h"
#include "Headers/Cell.h"
#include "Headers/Grid.h"
#include "Headers/Face.h"
#include "Headers/GravMass.h"
#include "Headers/onestep.h"
#include "Headers/IO.h"
#include "Headers/header.h"

int main(int argc, char **argv) {

  // start MPI 
  int error = mpiSetup(argc,argv);
  if (error==1) return(0);
  
  //grid
  struct Grid * theGrid = grid_create();
  grid_set_N_p(theGrid);
  grid_set_rz(theGrid);
  grid_set_Ncells_and_offset(theGrid);

  // gravMass
  struct GravMass *theGravMasses = gravMass_create(NP);
  gravMass_initialize(theGravMasses);
  gravMass_clean_pi(theGravMasses);
  
  // allocate memory for data 
  struct Cell ***theCells = cell_create(theGrid);
  cell_clean_pi(theCells,theGrid);
  // set initial data 
  cell_init(theCells,theGrid);
 
  //inter-processor syncs
  cell_syncproc_r(theCells,theGrid);
  cell_syncproc_z(theCells,theGrid);

  //set conserved quantities
  cell_prim2cons(theCells,theGrid);

  double dtcheck = T_MAX/NUM_CHECKPOINTS;
  double tcheck = dtcheck;
  
  int nfile=0;
  char filename[256];

  double t  = 0.0;
  while( t < T_MAX ){
    double dt = cell_mindt( theCells ,theGrid);
    printf("t: %e, dt: %e\n",t,dt);
    if( t+dt > T_MAX ) dt = T_MAX-t;
    cell_copy(theCells,theGrid);
    gravMass_copy(theGravMasses);
    onestep(theCells,theGrid,theGravMasses,dt,0.0);
    onestep(theCells,theGrid,theGravMasses,.5*dt,0.5);
    onestep_Psi(theCells,theGrid,dt);
    t += dt; 
    if( t>tcheck){
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
  cell_syncproc_r(theCells,theGrid);
  cell_syncproc_z(theCells,theGrid);

  // clean up
  cell_destroy(theCells,theGrid);
  grid_destroy(theGrid);
  gravMass_destroy(theGravMasses);

  // exit MPI 
  MPI_Barrier(grid_comm);
  MPI_Finalize();

  return(0);
}
