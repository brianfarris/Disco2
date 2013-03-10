#ifndef SIM_H
#define SIM_H
struct Sim;
struct Cell;
struct MPIsetup;

#ifdef SIM_PRIVATE_DEFS
struct Sim {
  double *r_faces;
  double *z_faces;
  int *N_p;
  int N_r_0;
  int N_z_0;
  int N_r_noghost;
  int N_z_noghost;
  int Nghost_rmin;
  int Nghost_rmax;
  int Nghost_zmin;
  int Nghost_zmax;
  int Ncells;
  int Ncells_global;
  int offset;
  //fixed parameters
  int InitialDataType;
  int GravMassType;
  int BoundTypeR;
  int BoundTypeZ;
  int Restart;
  int N_r_global;
  int N_z_global;
  int NP_CONST;
  double aspect;
  int ng;
  double T_MAX;
  int NUM_CHECKPOINTS;
  int NUM_DIAG_DUMP;
  int NUM_DIAG_MEASURE;
  double RMIN;
  double RMAX;
  double ZMIN;
  double ZMAX;
  int NUM_Q;
  int MOVE_CELLS;   
  int NumGravMass;
  double GAMMALAW;
  int INCLUDE_VISCOSITY;
  double EXPLICIT_VISCOSITY;
  double DIVB_CH;
  double DIVB_L;
  double CFL;
  double PLM;
  int POWELL;
  int GRAV2D;
  double G_EPS;
  double PHI_ORDER;
  double RHO_FLOOR;
  double CS_FLOOR;
  double CS_CAP;
  double VEL_CAP;
  int runtype;
  void (*single_init_ptr)(struct Cell ***,struct Sim *,int,int,int);
  void (*init_ptr)(struct Cell ***,struct Sim *);
};
#endif

//create and destroy
struct Sim *sim_create(struct MPIsetup * );
void sim_alloc_arr(struct Sim * , struct MPIsetup * );
void sim_destroy(struct Sim *); 
//access Sim data
double sim_RMIN(struct Sim * );
double sim_RMAX(struct Sim * );
double sim_ZMIN(struct Sim * );
double sim_ZMAX(struct Sim * );
int sim_N_r_0(struct Sim * );
int sim_N_z_0(struct Sim * );
int sim_N_p(struct Sim *,int);
double sim_r_faces(struct Sim *,int);
double sim_z_faces(struct Sim *,int);
int sim_N_r(struct Sim *);
int sim_N_z(struct Sim *);
int sim_Restart(struct Sim *);
int sim_BoundTypeR(struct Sim *);
int sim_BoundTypeZ(struct Sim *);
int sim_N_r_global(struct Sim *);
int sim_N_z_global(struct Sim *);
int sim_Ncells(struct Sim *);
int sim_Ncells_global(struct Sim *);
int sim_offset(struct Sim *);
int sim_ng(struct Sim *);
int sim_Nghost_rmin(struct Sim *);
int sim_Nghost_rmax(struct Sim *);
int sim_Nghost_zmin(struct Sim *);
int sim_Nghost_zmax(struct Sim *);
int sim_MOVE_CELLS(struct Sim *);
int sim_NumGravMass(struct Sim *);
double sim_GAMMALAW(struct Sim *);
int sim_INCLUDE_VISCOSITY(struct Sim *);
double sim_EXPLICIT_VISCOSITY(struct Sim *);
double sim_DIVB_CH(struct Sim *);
double sim_DIVB_L(struct Sim *);
double sim_CFL(struct Sim *);
double sim_PLM(struct Sim *);
int sim_POWELL(struct Sim *);
int sim_GRAV2D(struct Sim *);
double sim_G_EPS(struct Sim *);
double sim_PHI_ORDER(struct Sim *);
double sim_RHO_FLOOR(struct Sim *);
double sim_CS_FLOOR(struct Sim *);
double sim_CS_CAP(struct Sim *);
double sim_VEL_CAP(struct Sim *);
int sim_NUM_Q(struct Sim *);
double sim_get_T_MAX(struct Sim * );
int sim_NUM_CHECKPOINTS(struct Sim * );
int sim_NUM_DIAG_DUMP(struct Sim * );
int sim_NUM_DIAG_MEASURE(struct Sim * );
int sim_runtype(struct Sim * );
int sim_InitialDataType(struct Sim * );
int sim_GravMassType(struct Sim * );

//set Grid data
int sim_read_par_file(struct Sim * ,struct MPIsetup *, char * );
void sim_set_N_p(struct Sim *);
void sim_set_rz(struct Sim *,struct MPIsetup *);
void sim_set_misc(struct Sim *,struct MPIsetup *);
#endif
