#ifndef CELL_H
#define CELL_H
struct Cell;
struct Grid;
struct Face;
struct GravMass;
struct MPIsetup;

#ifdef CELL_PRIVATE_DEFS
struct Cell{
  double * prim;
  double * cons;
  double * RKcons;
  double * grad;
  double * gradp;
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
double *cell_prims(struct Cell *);
double *cell_grad(struct Cell *);
double *cell_gradp(struct Cell *);
struct Cell *cell_single(struct Cell ***,int,int,int);
double cell_tiph(struct Cell *);
double cell_dphi(struct Cell *);
double cell_wiph(struct Cell *);
//modify cell data
void cell_add_cons(struct Cell *, int, double);
void cell_add_divB(struct Cell *, double);
void cell_add_GradPsi(struct Cell *, int, double);
void cell_add_wiph(struct Cell *, double);
void cell_add_src( struct Cell *** ,struct Grid * , struct GravMass * , double );
void cell_add_visc_src( struct Cell *** ,struct Grid * , double );
void cell_mult_psi(struct Cell *, double);
void cell_clean_pi(struct Cell *** ,struct Grid *);
void cell_update_phi( struct Cell *** , struct Grid * , double , double );
void cell_update_dphi( struct Cell *** ,struct Grid * );
//clear
void cell_clear_w(struct Cell ***,struct Grid * );
void cell_clear_divB( struct Cell ***,struct Grid * );
void cell_clear_GradPsi( struct Cell ***,struct Grid * );
//set w
void cell_set_wcell(struct Cell ***,struct Grid *);
void cell_set_wrigid(struct Cell ***,struct Grid *);
//processor syncs
void cell_syncproc_r( struct Cell *** , struct Grid *,struct MPIsetup *);
void cell_syncproc_z( struct Cell *** , struct Grid *,struct MPIsetup *);
//PLM
void cell_plm_rz( struct Cell *** ,struct Grid *, struct Face * , int , int );
void cell_plm_p( struct Cell *** ,struct Grid * );
//boundary conditions
void cell_boundary_outflow_r( struct Cell *** , struct Face * ,struct Grid * ,struct MPIsetup *, int * );
void cell_boundary_outflow_z( struct Cell *** , struct Face * , struct Grid * ,struct MPIsetup *,int * );
void cell_boundary_fixed_r( struct Cell ***, struct Grid *,struct MPIsetup *);
//primitive-conservative conversion routines
void cell_calc_prim( struct Cell ***,struct Grid *);
void cell_prim2cons( double * , double * , double , double ,double,int );
void cell_calc_cons(struct Cell ***,struct Grid *); 
void cell_cons2prim( double * , double * , double , double ,struct Grid *,int );
//miscellaneous
double cell_mindt( struct Cell *** , struct Grid * );
void cell_copy(struct Cell ***,struct Grid * );
void cell_adjust_RK_cons( struct Cell *** , struct Grid * , double );
void cell_set_prim(struct Cell ***,int,int,int,int,double);//this will morph into a checkpoint restart routine
void cell_set_tiph(struct Cell ***,int,int,int,double);//this will morph into a checkpoint restart routine
void cell_print_cons(struct Cell *** ,struct Grid *,int);
void cell_print_prim(struct Cell ***  ,struct Grid *,int);
#endif
