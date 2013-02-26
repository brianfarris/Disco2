enum{RHO,PPP,URR,UPP,UZZ,BRR,BPP,BZZ,PSI};
enum{DDD,TAU,SRR,LLL,SZZ};
enum{C_FIXED,C_WCELL,C_WRIEMANN,C_RIGID};

#include "mpi.h"
//#define N_z_global  32 
//#define N_r_global  32
//#define ng          2
#define NUM_Q   9


//#define RMIN    1.0
//#define RMAX    4.0
//#define ZMIN   -1.0
//#define ZMAX    1.0
//#define T_MAX   50.0
//#define NUM_CHECKPOINTS 10
#define NP      1
#define Z_PERIODIC 1
#define GAMMALAW 1.66666666
#define INCLUDE_VISCOSITY 0
#define EXPLICIT_VISCOSITY 0.1
#define DIVB_CH 0.06
#define DIVB_L 0.1
#define CFL     0.5
#define PLM     1.0
#define MOVE_CELLS  C_RIGID
#define POWELL 1
#define GRAV2D 1
#define G_EPS 0.0
#define PHI_ORDER   2.0
#define RHO_FLOOR 0.00001
#define CS_FLOOR 0.0001
#define CS_CAP 1.0
#define VEL_CAP 10.0

/*
int MyProc,NumProcs;
int dim_MyProc[2];
int dim_NumProcs[2];
int left_Proc[2];
int right_Proc[2];
*/

MPI_Comm grid_comm;


