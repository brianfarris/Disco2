#include <stdio.h>
#include <stdlib.h>
#include "../Headers/Cell.h"
#include "../Headers/Sim.h"
#include "../Headers/header.h"
void (*cell_init_ptr(struct Sim * theSim))(struct Cell *** , struct Sim *,struct MPIsetup * ){
  if (sim_InitialDataType(theSim)==FLOCK){
    return(&cell_init_flock);
  } else if (sim_InitialDataType(theSim)==SHEAR){
    return(&cell_init_shear);
  } else if (sim_InitialDataType(theSim)==VORTEX){
    return(&cell_init_vortex);
  } else if (sim_InitialDataType(theSim)==STONE){
    return(&cell_init_stone);
  } else if (sim_InitialDataType(theSim)==FIELDLOOP){
    return(&cell_init_fieldloop);
  } else if (sim_InitialDataType(theSim)==PSIGRAD){
    return(&cell_init_psigrad);
  } else if (sim_InitialDataType(theSim)==TORUS){
    return(&cell_init_torus);
  } else if (sim_InitialDataType(theSim)==MILOS_MACFADYEN){
    return(&cell_init_milos_macfadyen);
  } else if (sim_InitialDataType(theSim)==MHDEXP){
    return(&cell_init_mhdexp);
  } else if (sim_InitialDataType(theSim)==VISCRING){
    return(&cell_init_viscring_cstnu);
  } else if (sim_InitialDataType(theSim)==TESTING){
    return(&cell_init_testing);
  } else if (sim_InitialDataType(theSim)==DIVERGENCE){ 
    return(&cell_init_divergence);
  } else if (sim_InitialDataType(theSim)==RAD_DOM){
    return(&cell_init_rad_dom);
  } else if (sim_InitialDataType(theSim)==MIDDLE){
    return(&cell_init_middle);
  } else{
    printf("ERROR, not a valid global initial data routine choice\n");
    exit(0);
  }
}

void (*cell_single_init_ptr(struct Sim * theSim))(struct Cell * , struct Sim *,int,int,int ){
  if (sim_InitialDataType(theSim)==FLOCK){
    return(&cell_single_init_flock);
  } else if (sim_InitialDataType(theSim)==SHEAR){
    return(&cell_single_init_shear);
  } else if (sim_InitialDataType(theSim)==VORTEX){
    return(&cell_single_init_vortex);
  } else if (sim_InitialDataType(theSim)==STONE){
    return(&cell_single_init_stone);
  } else if (sim_InitialDataType(theSim)==FIELDLOOP){
    return(&cell_single_init_fieldloop);
  } else if (sim_InitialDataType(theSim)==PSIGRAD){
    return(&cell_single_init_psigrad);
  } else if (sim_InitialDataType(theSim)==TORUS){
    return(&cell_single_init_torus);
  } else if (sim_InitialDataType(theSim)==MILOS_MACFADYEN){
    return(&cell_single_init_milos_macfadyen);
  } else if (sim_InitialDataType(theSim)==MHDEXP){
    return(&cell_single_init_mhdexp);
  } else if (sim_InitialDataType(theSim)==VISCRING){
    return(&cell_single_init_viscring_cstnu);
  } else if (sim_InitialDataType(theSim)==TESTING){
    return(&cell_single_init_testing);
  } else if (sim_InitialDataType(theSim)==DIVERGENCE){
    return(&cell_single_init_divergence);
  } else if (sim_InitialDataType(theSim)==MIDDLE){
    return(&cell_single_init_middle);
  } else{
    printf("ERROR, not a valid single point initial data routine choice\n");
    exit(0);
  }
}


