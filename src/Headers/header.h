enum{RHO,PPP,URR,UPP,UZZ,BRR,BPP,BZZ,PSI,ARR,APP,AZZ,PHI};
enum{DDD,TAU,SRR,LLL,SZZ};
enum{C_FIXED,C_WCELL,C_WRIEMANN,C_RIGID,C_KEPLER,C_OMEGA1,C_MILOS};
enum{A_NONE,A_OMEGA20,A_KEPLER,A_WEIRD};
enum{LEFT,LEFTSTAR,RIGHTSTAR,RIGHT};
enum{EULER,MHD};
enum{VAR_INT,VAR_DOUB,VAR_STR};
enum{BOUND_FIXED,BOUND_OUTFLOW,BOUND_PERIODIC};
enum{FLOCK,SHEAR,VORTEX,STONE,FIELDLOOP,PSIGRAD,TORUS,MILOS_MACFADYEN,MHDEXP,VISCRING,TESTING};
enum{NONE,SINGLE,BINARY};
enum{LINEAR,SPLINE};
//unify the 2 below at some point
enum{R_DIR,Z_DIR};
enum{RDIRECTION,PDIRECTION,ZDIRECTION};
enum{HLLC,HLL};
#include "mpi.h"
MPI_Comm sim_comm;
#define farris_mpi_factorization 0
#define NO_W_IN_CFL 1
#define KEP_BNDRY 0
#define BzZ 0
#define BNORM_AVG 0
#define VISC_CONST 0
#define VISC_OLD 0
#define INCLUDE_ALL_VISC_TERMS 0
#define CHECKPOINTING 
#define TVISC_FAC 1.0
double time_global;
