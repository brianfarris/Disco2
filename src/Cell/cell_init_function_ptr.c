#include <stdio.h>
#include <stdlib.h>
#include "../Headers/Cell.h"
#include "../Headers/Sim.h"
#include "../Headers/header.h"
void (*cell_init_ptr(struct Sim * theSim))(struct Cell *** , struct Sim *,struct GravMass * ,struct MPIsetup * ){
    if (sim_InitialDataType(theSim)==SHEAR){
        return(&cell_init_shear);
    } else if (sim_InitialDataType(theSim)==VORTEX){
        return(&cell_init_vortex);
    } else if (sim_InitialDataType(theSim)==TORUS){
        return(&cell_init_torus);
    } else if (sim_InitialDataType(theSim)==MILOS_MACFADYEN){
        return(&cell_init_milos_macfadyen);
    } else if (sim_InitialDataType(theSim)==RAD_DOM){
        return(&cell_init_rad_dom);
    } else if (sim_InitialDataType(theSim)==MIDDLE){
        return(&cell_init_middle);
    } else if (sim_InitialDataType(theSim)==SSTEST){
        return(&cell_init_SStest);
    } else if (sim_InitialDataType(theSim)==BETA){
        return(&cell_init_beta);
     } else{
        printf("ERROR, not a valid global initial data routine choice\n");
        exit(0);
    }
}

void (*cell_single_init_ptr(struct Sim * theSim))(struct Cell * , struct Sim *,struct GravMass * ,int,int,int ){
    if (sim_InitialDataType(theSim)==SHEAR){
        return(&cell_single_init_shear);
    } else if (sim_InitialDataType(theSim)==VORTEX){
        return(&cell_single_init_vortex);
    } else if (sim_InitialDataType(theSim)==TORUS){
        return(&cell_single_init_torus);
    } else if (sim_InitialDataType(theSim)==MILOS_MACFADYEN){
        return(&cell_single_init_milos_macfadyen);
    } else if (sim_InitialDataType(theSim)==RAD_DOM){
        return(&cell_single_init_rad_dom);
     } else if (sim_InitialDataType(theSim)==MIDDLE){
        return(&cell_single_init_middle);
    } else if (sim_InitialDataType(theSim)==SSTEST){
        return(&cell_single_init_SStest);
    } else if (sim_InitialDataType(theSim)==BETA){
        return(&cell_single_init_beta);
     } else{
        printf("ERROR, not a valid single point initial data routine choice\n");
        exit(0);
    }
}


