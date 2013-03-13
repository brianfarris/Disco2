#define TIMESTEP_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include "../Headers/Sim.h"
#include "../Headers/MPIsetup.h"
#include "../Headers/GravMass.h"
#include "../Headers/Cell.h"
#include "../Headers/TimeStep.h"
#include "../Headers/header.h"

void timestep_rk2(struct TimeStep * theTimeStep, struct Sim * theSim,
    struct Cell *** theCells, struct GravMass * theGravMasses, struct MPIsetup * theMPIsetup) {

  cell_clear_w(theCells,theSim);
  cell_set( theCells ,theSim);
  timestep_set_dt(theTimeStep,theCells,theSim);
  cell_copy(theCells,theSim);
  gravMass_copy(theGravMasses,theSim);
  timestep_set_RK(theTimeStep,0.0);
  timestep_substep(theTimeStep,theCells,theSim,theGravMasses,theMPIsetup,1.0);
  gravMass_move(theGravMasses,1.0*timestep_dt(theTimeStep));
  timestep_set_RK(theTimeStep,0.5);
  timestep_substep(theTimeStep,theCells,theSim,theGravMasses,theMPIsetup,0.5);
  gravMass_move(theGravMasses,0.5*timestep_dt(theTimeStep));
  if (sim_runtype(theSim)==1) timestep_update_Psi(theTimeStep,theCells,theSim,theMPIsetup);
  timestep_update_t(theTimeStep); 

}

void timestep_forward_euler(struct TimeStep * theTimeStep, struct Sim * theSim,
    struct Cell *** theCells, struct GravMass * theGravMasses, struct MPIsetup * theMPIsetup) {

  cell_clear_w(theCells,theSim);
  cell_set( theCells ,theSim);
  timestep_set_dt(theTimeStep,theCells,theSim);
  cell_copy(theCells,theSim);
  gravMass_copy(theGravMasses,theSim);
  timestep_set_RK(theTimeStep,0.0);
  timestep_substep(theTimeStep,theCells,theSim,theGravMasses,theMPIsetup,1.0);
  gravMass_move(theGravMasses,1.0*timestep_dt(theTimeStep));
  if (sim_runtype(theSim)==1) timestep_update_Psi(theTimeStep,theCells,theSim,theMPIsetup);
  timestep_update_t(theTimeStep); 
  
}
