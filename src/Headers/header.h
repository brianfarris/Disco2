enum{RHO,PPP,URR,UPP,UZZ,BRR,BPP,BZZ,PSI};
enum{DDD,TAU,SRR,LLL,SZZ};
enum{C_FIXED,C_WCELL,C_WRIEMANN,C_RIGID,C_KEPLER,C_OMEGA20,C_MILOS};
enum{A_NONE,A_OMEGA20,A_KEPLER,A_WEIRD};
enum{LEFT,LEFTSTAR,RIGHTSTAR,RIGHT};
enum{EULER,MHD};
enum{VAR_INT,VAR_DOUB,VAR_STR};
enum{BOUND_FIXED,BOUND_OUTFLOW,BOUND_PERIODIC};
enum{FLOCK,SHEAR,VORTEX,STONE,FIELDLOOP,PSIGRAD,TORUS,MILOS_MACFADYEN,MHDEXP,VISCRING, TYPEIISD_MIGRATION, COOLTEST};
enum{NONE,SINGLE,BINARY,LIVEBINARY};
//unify the 2 below at some point
enum{R_DIR,Z_DIR};
enum{RDIRECTION,PDIRECTION,ZDIRECTION};
enum{HLLC,HLL};
#include "mpi.h"
MPI_Comm sim_comm;
#define farris_mpi_factorization 0
#define NO_W_IN_CFL 0  // if 1 then no w in cfl
#define KEP_BNDRY 0
#define BzZ 0
#define BNORM_AVG 1
#define VISC_CONST 1   //if ==1 then nu=cst else alpha law unless IC=shear - see Riemann_routines.c
#define VISC_OLD 0     // 0 off 1 on
#define INCLUDE_ALL_VISC_TERMS 1  // 0 -> only rphi terms 1 all terms
//#define CHECKPOINTING //Uncomment to use checkpointing in hdf5 parallel installed!
double time_global;
