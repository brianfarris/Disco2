#define MPISETUP_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include "../Headers/MPIsetup.h"
#include "../Headers/header.h"

int mpisetup_check_rin_bndry(struct MPIsetup * theMPIsetup){
  return(theMPIsetup->dim_MyProc[0]==0);
}
int mpisetup_check_rout_bndry(struct MPIsetup * theMPIsetup){
  return(theMPIsetup->dim_MyProc[0]== (theMPIsetup->dim_NumProcs[0]-1));
}
int mpisetup_check_zbot_bndry(struct MPIsetup * theMPIsetup){
  return(theMPIsetup->dim_MyProc[1]==1);
}
int mpisetup_check_ztop_bndry(struct MPIsetup * theMPIsetup){
  return(theMPIsetup->dim_MyProc[1]== (theMPIsetup->dim_NumProcs[1]-1));
}
int mpisetup_MyProc(struct MPIsetup * theMPIsetup){
  return(theMPIsetup->MyProc);
}
int mpisetup_NumProcs(struct MPIsetup * theMPIsetup){
  return(theMPIsetup->NumProcs);
}
int mpisetup_dim_MyProc(struct MPIsetup * theMPIsetup,int dim){
  return(theMPIsetup->dim_MyProc[dim]);
}
int mpisetup_dim_NumProcs(struct MPIsetup * theMPIsetup,int direction){
  return(theMPIsetup->dim_NumProcs[direction]);
}
int mpisetup_left_Proc(struct MPIsetup * theMPIsetup,int direction){
  return(theMPIsetup->left_Proc[direction]);
}
int mpisetup_right_Proc(struct MPIsetup * theMPIsetup,int direction){
  return(theMPIsetup->right_Proc[direction]);
}

