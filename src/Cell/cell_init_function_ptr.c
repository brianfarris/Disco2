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
  } else if (sim_InitialDataType(theSim)==UNIFORM){
    return(&cell_init_uniform);
  } else if (sim_InitialDataType(theSim)==SHOCK1){
    return(&cell_init_shock1);
  } else if (sim_InitialDataType(theSim)==SHOCK2){
    return(&cell_init_shock2);
  } else if (sim_InitialDataType(theSim)==SHOCK3){
    return(&cell_init_shock3);
  } else if (sim_InitialDataType(theSim)==SHOCK4){
    return(&cell_init_shock4);
  } else if (sim_InitialDataType(theSim)==ISENTROPE){
    return(&cell_init_isentrope);
  } else if (sim_InitialDataType(theSim)==GBONDI){
    return(&cell_init_gbondi);
  } else if (sim_InitialDataType(theSim)==GBONDI2){
    return(&cell_init_gbondi2);
  } else if (sim_InitialDataType(theSim)==EQUIL1){
    return(&cell_init_equil1);
  } else if (sim_InitialDataType(theSim)==EQUIL2){
    return(&cell_init_equil2);
  } else if (sim_InitialDataType(theSim)==SSDISC){
    return(&cell_init_ssdisc);
  } else if (sim_InitialDataType(theSim)==NTDISC){
    return(&cell_init_ntdisc);
  } else if (sim_InitialDataType(theSim)==CNSTDISC){
    return(&cell_init_cnstdisc);
  } else if (sim_InitialDataType(theSim)==CARTSHEAR){
    return(&cell_init_cartshear);
  } else if (sim_InitialDataType(theSim)==DISCTEST){
    return(&cell_init_disctest);
  } else if (sim_InitialDataType(theSim)==ACCDISC){
    return(&cell_init_accdisc);
  } else if (sim_InitialDataType(theSim)==ADAF){
    return(&cell_init_adaf);
  } else if (sim_InitialDataType(theSim)==ADAF2){
    return(&cell_init_adaf2);
  } else if (sim_InitialDataType(theSim)==ATMO){
    return(&cell_init_atmo);
  } else if (sim_InitialDataType(theSim)==ORBIT){
    return(&cell_init_orbit);
  } else if (sim_InitialDataType(theSim)==MINIDISC){
    return(&cell_init_minidisc);
  } else if (sim_InitialDataType(theSim)==KICK){
    return(&cell_init_kick);
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
   } else if (sim_InitialDataType(theSim)==UNIFORM){
    return(&cell_single_init_uniform);
   } else if (sim_InitialDataType(theSim)==SHOCK1){
    return(&cell_single_init_shock1);
   } else if (sim_InitialDataType(theSim)==SHOCK2){
    return(&cell_single_init_shock2);
   } else if (sim_InitialDataType(theSim)==SHOCK3){
    return(&cell_single_init_shock3);
   } else if (sim_InitialDataType(theSim)==SHOCK4){
    return(&cell_single_init_shock4);
   } else if (sim_InitialDataType(theSim)==ISENTROPE){
    return(&cell_single_init_isentrope);
   } else if (sim_InitialDataType(theSim)==GBONDI){
    return(&cell_single_init_gbondi);
   } else if (sim_InitialDataType(theSim)==GBONDI2){
    return(&cell_single_init_gbondi2);
   } else if (sim_InitialDataType(theSim)==EQUIL1){
    return(&cell_single_init_equil1);
   } else if (sim_InitialDataType(theSim)==EQUIL2){
    return(&cell_single_init_equil2);
   } else if (sim_InitialDataType(theSim)==SSDISC){
    return(&cell_single_init_ssdisc);
   } else if (sim_InitialDataType(theSim)==NTDISC){
    return(&cell_single_init_ntdisc);
   } else if (sim_InitialDataType(theSim)==CNSTDISC){
    return(&cell_single_init_cnstdisc);
   } else if (sim_InitialDataType(theSim)==CARTSHEAR){
    return(&cell_single_init_cartshear);
   } else if (sim_InitialDataType(theSim)==DISCTEST){
    return(&cell_single_init_disctest);
   } else if (sim_InitialDataType(theSim)==ACCDISC){
    return(&cell_single_init_accdisc);
   } else if (sim_InitialDataType(theSim)==ADAF){
    return(&cell_single_init_adaf);
   } else if (sim_InitialDataType(theSim)==ADAF2){
    return(&cell_single_init_adaf2);
   } else if (sim_InitialDataType(theSim)==ATMO){
    return(&cell_single_init_atmo);
   } else if (sim_InitialDataType(theSim)==ORBIT){
    return(&cell_single_init_orbit);
   } else if (sim_InitialDataType(theSim)==MINIDISC){
    return(&cell_single_init_minidisc);
   } else if (sim_InitialDataType(theSim)==KICK){
    return(&cell_single_init_kick);
  } else{
    printf("ERROR: Do not recognize initial data selection.\n");
    exit(0);
  }
}


