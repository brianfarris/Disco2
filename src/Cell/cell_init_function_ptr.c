#include <stdio.h>
#include <stdlib.h>
#include "../Headers/Cell.h"
#include "../Headers/Sim.h"
#include "../Headers/header.h"
void (*cell_init_ptr(struct Sim * theSim))(struct Cell *** , struct Sim *,struct MPIsetup * ){
  if (sim_InitialDataType(theSim)==VORTEX){
    return(&cell_init_vortex);
  } else if (sim_InitialDataType(theSim)==TORUS){
    return(&cell_init_torus);
  } else if (sim_InitialDataType(theSim)==BONDI){
    return(&cell_init_bondi);
   } else{
    printf("ERROR\n");
    exit(0);
  }
}

void (*cell_single_init_ptr(struct Sim * theSim))(struct Cell * , struct Sim *,int,int,int ){
  if (sim_InitialDataType(theSim)==VORTEX){
    return(&cell_single_init_vortex);
  } else if (sim_InitialDataType(theSim)==TORUS){
    return(&cell_single_init_torus);
   } else if (sim_InitialDataType(theSim)==BONDI){
    return(&cell_single_init_bondi);
  } else{
    printf("ERROR\n");
    exit(0);
  }
}


