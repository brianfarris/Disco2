#define IO_PRIVATE_DEFS
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "../Headers/IO.h"
#include "../Headers/TimeStep.h"
#include "../Headers/Sim.h"
#include "../Headers/header.h"

void io_setup(struct IO * theIO,struct Sim * theSim,struct TimeStep * theTimeStep){
  theIO->nfile=0;
  if(sim_NUM_CHECKPOINTS(theSim) < 1)
    theIO->dtcheck = sim_get_T_MAX(theSim);
  else
    theIO->dtcheck = sim_get_T_MAX(theSim)/(sim_NUM_CHECKPOINTS(theSim)-1);
  theIO->tcheck = 0.0;
  while( theIO->tcheck < timestep_get_t(theTimeStep)){
    theIO->tcheck += theIO->dtcheck;
    ++theIO->nfile;
  }
  sprintf(theIO->filename,"checkpoint_%04d.h5",io_nfile(theIO));
}

