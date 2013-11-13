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
  } else if (sim_InitialDataType(theSim)==SHOCK1){
    return(&cell_init_shock1);
  } else if (sim_InitialDataType(theSim)==UNIFORM){
    return(&cell_init_uniform);
   } else{
    printf("ERROR: Do not recognize initial data selection.\n");
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
   } else if (sim_InitialDataType(theSim)==SHOCK1){
    return(&cell_single_init_shock1);
   } else if (sim_InitialDataType(theSim)==UNIFORM){
    return(&cell_single_init_uniform);
  } else{
    printf("ERROR: Do not recognize initial data selection.\n");
    exit(0);
  }
}


