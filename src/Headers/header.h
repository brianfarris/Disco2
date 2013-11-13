enum{RHO,PPP,URR,UPP,UZZ,BRR,BPP,BZZ,PSI};
enum{DDD,TAU,SRR,LLL,SZZ};
enum{C_FIXED,C_WCELL,C_WRIEMANN,C_RIGID,C_KEPLER,C_OMEGA20,C_MILOS};
enum{A_NONE,A_OMEGA20,A_KEPLER,A_WEIRD};
enum{LEFT,LEFTSTAR,RIGHTSTAR,RIGHT};
enum{VAR_INT,VAR_DOUB,VAR_STR};
enum{BOUND_FIXED,BOUND_OUTFLOW,BOUND_PERIODIC};
enum{VORTEX,TORUS,BONDI,SHOCK1,UNIFORM};
enum{NONE,SINGLE,BINARY};
//unify the 2 below at some point
enum{R_DIR,Z_DIR};
enum{RDIRECTION,PDIRECTION,ZDIRECTION};
enum{HLLC,HLL};
//Background Types
enum{NEWTON,GR};
//Metric Types
enum{SR,SCHWARZSCHILD,KERR};
#include "mpi.h"
MPI_Comm sim_comm;
#define farris_mpi_factorization 0
#define KEP_BNDRY 0
//#define NPCAP 256 
double time_global;
