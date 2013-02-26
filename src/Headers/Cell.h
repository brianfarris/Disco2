#ifndef CELL_H
#define CELL_H
#define NUM_Q   9
struct Cell;
struct Grid;
struct Face;
struct GravMass;
struct MPIsetup;

#ifdef CELL_PRIVATE_DEFS
struct Cell{
  double prim[NUM_Q];
  double cons[NUM_Q];
  double RKcons[NUM_Q];
  double grad[NUM_Q];
  double gradp[NUM_Q];
  double tiph;
  double RKtiph;
  double dphi;
  double wiph;
  double divB;
  double GradPsi[3];
};
#endif
//create and destroy
struct Cell ***cell_create(struct Grid *,struct MPIsetup *);
void cell_destroy(struct Cell ***,struct Grid *);
//initial data
void cell_init(struct Cell ***,struct Grid *,struct MPIsetup *);
void cell_single_init(struct Cell ***, struct Grid *,int ,int ,int );
//retrieve data
double *cell_get_prims(struct Cell ***,int ,int ,int );//redundant
double *cell_single_get_prims(struct Cell *);//redundant
double *cell_get_cons(struct Cell ***,int ,int ,int );
double *cell_single_get_grad(struct Cell *);
double *cell_single_get_gradp(struct Cell *);
struct Cell *cell_pointer(struct Cell ***,int,int,int);
double cell_tiph(struct Cell ***,int,int,int);//redundant
double cell_single_tiph(struct Cell *);//redundant
double cell_dphi(struct Cell ***,int,int,int);//redundant
double cell_single_dphi(struct Cell *);//redundant
double cell_single_wiph(struct Cell *);
//modify cell data
void cell_add_cons(struct Cell *, int, double);
void cell_add_divB(struct Cell *, double);
void cell_add_GradPsi(struct Cell *, int, double);
void cell_add_wiph(struct Cell *, double);
void cell_add_src( struct Cell *** ,struct Grid * , struct GravMass * , double );
void cell_mult_psi(struct Cell *, double);
void cell_clean_pi(struct Cell *** ,struct Grid *);
void cell_set_tiph(struct Cell ***,int,int,int,double);
void cell_set_wcell(struct Cell ***,struct Grid *);
void cell_set_wrigid(struct Cell ***,struct Grid *);
void cell_clear_w(struct Cell ***,struct Grid * );
void cell_clear_divB( struct Cell ***,struct Grid * );
void cell_clear_GradPsi( struct Cell ***,struct Grid * );
void cell_update_phi( struct Cell *** , struct Grid * , double , double );
void cell_update_dphi( struct Cell *** ,struct Grid * );
//processor syncs
void cell_syncproc_r( struct Cell *** , struct Grid *,struct MPIsetup *);
void cell_syncproc_z( struct Cell *** , struct Grid *,struct MPIsetup *);
//riemann
void cell_riemann_p( struct Cell * , struct Cell * , double , double , double );
void cell_plm_rz( struct Cell *** ,struct Grid *, struct Face * , int , int );
void cell_plm_p( struct Cell *** ,struct Grid * );
void cell_flux_p( struct Cell *** ,struct Grid *, double );
//boundary conditions
void cell_boundary_outflow_r( struct Cell *** , struct Face * ,struct Grid * ,struct MPIsetup *, int * );
void cell_boundary_outflow_z( struct Cell *** , struct Face * , struct Grid * ,struct MPIsetup *,int * );
void cell_boundary_fixed_r( struct Cell ***, struct Grid *,struct MPIsetup *);
//primitive-conservative conversion routines
void cell_calc_prim( struct Cell ***,struct Grid *);
void cell_prim2cons(struct Cell ***,struct Grid *); 
//miscellaneous
double cell_mindt( struct Cell *** , struct Grid * );
void cell_copy(struct Cell ***,struct Grid * );
void cell_adjust_RK_cons( struct Cell *** , struct Grid * , double );
void cell_set_prim(struct Cell ***,int,int,int,int,double);//this will morph into a checkpoint restart routine

//void cell_printscreen(struct Cell ***,struct Grid * );
#endif
