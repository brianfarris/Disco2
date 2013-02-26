#ifndef MPISETUP_H
#define MPISETUP_H
struct MPIsetup;

#ifdef MPISETUP_PRIVATE_DEFS
struct MPIsetup {
  int ndims_mpi;
  int wraparound[2];
  int reorder;
  int MyProc;
  int NumProcs;
  int dim_MyProc[2];
  int dim_NumProcs[2];
  int left_Proc[2];
  int right_Proc[2];
};
#endif

//create and destroy
struct MPIsetup *mpisetup_create(int,char **);
void mpisetup_destroy(struct MPIsetup *); 
//set MPI info
void mpisetup_setprocs(struct MPIsetup * );
void mpisetup_cart_create(struct MPIsetup * );
void mpisetup_left_right(struct MPIsetup * );
//check whether we are on a boundary proc
int mpisetup_check_rin_bndry(struct MPIsetup * );
int mpisetup_check_rout_bndry(struct MPIsetup * );
int mpisetup_check_zbot_bndry(struct MPIsetup * );
int mpisetup_check_ztop_bndry(struct MPIsetup * );
//access MPI info
int mpisetup_MyProc(struct MPIsetup * );
int mpisetup_NumProcs(struct MPIsetup * );
int * mpisetup_dim_MyProc(struct MPIsetup * );
int * mpisetup_dim_NumProcs(struct MPIsetup * );
int * mpisetup_left_Proc(struct MPIsetup * );
int * mpisetup_right_Proc(struct MPIsetup * );
#endif
