#define SIM_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include "../Headers/Sim.h"
#include "../Headers/header.h"


double sim_MIN(struct Sim * theSim,int direction){
  return(theSim->MIN[direction]);
}
double sim_MAX(struct Sim * theSim,int direction){
  return(theSim->MAX[direction]);
}
int sim_N0(struct Sim * theSim, int direction){
  return(theSim->N0[direction]);
}
int sim_N_p(struct Sim *theSim,int i){
  return(theSim->N_p[i]);
}
//double sim_r_faces(struct Sim *theSim,int i){
//  return(theSim->r_faces[i+1]);
//}
//double sim_z_faces(struct Sim *theSim,int k){
//  return(theSim->z_faces[k+1]);
//}
double sim_FacePos(struct Sim *theSim,int index,int direction){
  if (direction==R_DIR){
    //printf("index: %d, theSim->r_faces[index+1]: %e\n",index,theSim->r_faces[index+1]);
    return(theSim->r_faces[index+1]);
  } else if (direction==Z_DIR){
    return(theSim->z_faces[index+1]);
  } else{
    printf("ERROR\n");
    exit(0);
  }
}

int sim_N(struct Sim *theSim, int direction){
  return(theSim->N_noghost[direction]+theSim->Nghost_min[direction]+theSim->Nghost_max[direction]);
}
int sim_BoundTypeR(struct Sim *theSim){
  return(theSim->BoundTypeR);
}
int sim_BoundTypeZ(struct Sim *theSim){
  return(theSim->BoundTypeZ);
}
int sim_Ncells(struct Sim *theSim){
  return(theSim->Ncells);
}
int sim_NoInnerBC(struct Sim *theSim){
  return(theSim->NoInnerBC);
}
int sim_Restart(struct Sim *theSim){
  return(theSim->Restart);
}
int sim_Ncells_global(struct Sim *theSim){
  return(theSim->Ncells_global);
}
int sim_offset(struct Sim *theSim){
  return(theSim->offset);
}
int sim_Nghost_min(struct Sim *theSim,int direction){
  return(theSim->Nghost_min[direction]);
}
int sim_Nghost_max(struct Sim *theSim,int direction){
  return(theSim->Nghost_max[direction]);
}
int sim_ng(struct Sim *theSim){
  return(theSim->ng);
}
int sim_N_global(struct Sim *theSim,int direction){
  return(theSim->N_global[direction]);
}
int sim_NUM_Q(struct Sim *theSim){
  return(theSim->NUM_Q);
}
int sim_NUM_C(struct Sim *theSim){
  return(theSim->NUM_C);
}
int sim_W_NUMERIC_TYPE(struct Sim *theSim){
  return(theSim->W_NUMERIC_TYPE);
}
int sim_NumGravMass(struct Sim *theSim){
  return(theSim->NumGravMass);
}
double sim_MassRatio(struct Sim *theSim){
  return(theSim->MassRatio);
}
double sim_OrbShrinkTscale(struct Sim *theSim){
  return(theSim->OrbShrinkTscale);
}
double sim_OrbShrinkT0(struct Sim *theSim){
  return(theSim->OrbShrinkT0);
}
int sim_Riemann(struct Sim *theSim){
  return(theSim->Riemann);
}
double sim_GAMMALAW(struct Sim *theSim){
  return(theSim->GAMMALAW);
}
int sim_COOLING(struct Sim *theSim){
  return(theSim->COOLING);
}
double sim_EXPLICIT_VISCOSITY(struct Sim *theSim){
  return(theSim->EXPLICIT_VISCOSITY);
}
int sim_VISC_CONST(struct Sim *theSim){
  return(theSim->VISC_CONST);
}
double sim_PoRho_r1(struct Sim *theSim){
    return(theSim->PoRho_r1);
}
double sim_CFL(struct Sim *theSim){
  return(theSim->CFL);
}
double sim_PLM(struct Sim *theSim){
  return(theSim->PLM);
}
int sim_W_ANALYTIC_TYPE(struct Sim *theSim){
  return(theSim->W_ANALYTIC_TYPE);
}
int sim_GRAV2D(struct Sim *theSim){
  return(theSim->GRAV2D);
}
double sim_G_EPS(struct Sim *theSim){
  return(theSim->G_EPS);
}
int sim_RhoSinkOn(struct Sim *theSim){
  return(theSim->RhoSinkOn);
}
int sim_PHI_ORDER(struct Sim *theSim){
  return(theSim->PHI_ORDER);
}
double sim_RHO_FLOOR(struct Sim *theSim) {
  return(theSim->RHO_FLOOR);
}
double sim_CS_FLOOR(struct Sim *theSim){
  return(theSim->CS_FLOOR);
}
double sim_CS_CAP(struct Sim *theSim){
  return(theSim->CS_CAP);
}
double sim_VEL_CAP(struct Sim *theSim){
  return(theSim->VEL_CAP);
}
double sim_get_T_MAX(struct Sim * theSim){
  return(theSim->T_MAX);
}
int sim_NUM_CHECKPOINTS(struct Sim * theSim){
  return(theSim->NUM_CHECKPOINTS);
}
int sim_NUM_DIAG_DUMP(struct Sim * theSim){
  return(theSim->NUM_DIAG_DUMP);
}
int sim_NUM_DIAG_MEASURE(struct Sim * theSim){
  return(theSim->NUM_DIAG_MEASURE);
}
int sim_InitialDataType(struct Sim * theSim){
  return(theSim->InitialDataType);
}
int sim_GravMassType(struct Sim * theSim){
  return(theSim->GravMassType);
}
double sim_DAMP_TIME(struct Sim * theSim){
  return(theSim->DAMP_TIME);
}
double sim_RDAMP_INNER(struct Sim * theSim){
  return(theSim->RDAMP_INNER);
}
double sim_RDAMP_OUTER(struct Sim * theSim){
  return(theSim->RDAMP_OUTER);
}

