#ifndef GRID_H
#define GRID_H
struct Grid;
struct MPIsetup;

#ifdef GRID_PRIVATE_DEFS
struct Grid {
  double *r_faces;
  double *z_faces;
  int *N_p;
  int N_r;
  int N_z;
  int Nghost_rmin;
  int Nghost_rmax;
  int Nghost_zmin;
  int Nghost_zmax;
  int Ncells;
  int Ncells_global;
  int offset;
  int Restart;
  int N_r_global;
  int N_z_global;
  int ng;
  double T_MAX;
  int NUM_CHECKPOINTS;
  double RMIN;
  double RMAX;
  double ZMIN;
  double ZMAX;
  int NUM_Q;
  int Z_PERIODIC;
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
};
#endif

//create and destroy
struct Grid *grid_create(struct MPIsetup * );
void grid_destroy(struct Grid *); 
//access grid data
int grid_N_p(struct Grid *,int);
double grid_r_faces(struct Grid *,int);
double grid_z_faces(struct Grid *,int);
int grid_N_r(struct Grid *);
int grid_N_z(struct Grid *);
int grid_Restart(struct Grid *);
int grid_N_z_global(struct Grid *);
int grid_Ncells(struct Grid *);
int grid_Ncells_global(struct Grid *);
int grid_offset(struct Grid *);
int grid_ng(struct Grid *);
int grid_Nghost_rmin(struct Grid *);
int grid_Nghost_rmax(struct Grid *);
int grid_Nghost_zmin(struct Grid *);
int grid_Nghost_zmax(struct Grid *);
int grid_MOVE_CELLS(struct Grid *);
int grid_NumGravMass(struct Grid *);
double grid_GAMMALAW(struct Grid *);
int grid_INCLUDE_VISCOSITY(struct Grid *);
double grid_EXPLICIT_VISCOSITY(struct Grid *);
double grid_DIVB_CH(struct Grid *);
double grid_DIVB_L(struct Grid *);
double grid_CFL(struct Grid *);
double grid_PLM(struct Grid *);
int grid_POWELL(struct Grid *);
int grid_GRAV2D(struct Grid *);
double grid_G_EPS(struct Grid *);
double grid_PHI_ORDER(struct Grid *);
double grid_RHO_FLOOR(struct Grid *);
double grid_CS_FLOOR(struct Grid *);
double grid_CS_CAP(struct Grid *);
double grid_VEL_CAP(struct Grid *);
int grid_NUM_Q(struct Grid *);
double grid_get_T_MAX(struct Grid * );
double grid_NUM_CHECKPOINTS(struct Grid * );
double grid_runtype(struct Grid * );


//set grid data
int grid_read_par_file(struct Grid * ,struct MPIsetup *, char * );
void grid_set_N_p(struct Grid *);
void grid_set_rz(struct Grid *,struct MPIsetup *);
void grid_set_Ncells_and_offset(struct Grid *,struct MPIsetup *);
#endif
