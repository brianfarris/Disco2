enum{RHO,PPP,URR,UPP,UZZ,BRR,BPP,BZZ,PSI};
enum{DDD,TAU,SRR,LLL,SZZ};
enum{C_FIXED,C_WCELL,C_WRIEMANN,C_RIGID};
enum{LEFT,LEFTSTAR,RIGHTSTAR,RIGHT};
enum{EULER,MHD};
enum{VAR_INT,VAR_DOUB,VAR_STR};
enum{BOUND_FIXED,BOUND_OUTFLOW,BOUND_PERIODIC};
enum{FLOCK,SHEAR};
enum{NONE,SINGLE,BINARY};
#include "mpi.h"
MPI_Comm grid_comm;


