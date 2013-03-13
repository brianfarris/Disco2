#ifndef CELL_H
#define CELL_H
struct Cell;
struct Sim;
struct Face;
struct GravMass;
struct MPIsetup;
struct TimeStep;

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
struct Cell ***cell_create(struct Sim *,struct MPIsetup *);
void cell_destroy(struct Cell ***,struct Sim *);
struct Cell * cell_single_create(struct Sim * );
void cell_single_destroy(struct Cell * );
//initial data
void cell_init_flock(struct Cell ***,struct Sim *, struct MPIsetup *);
void cell_single_init_flock(struct Cell *, struct Sim *,int ,int ,int );
void cell_init_shear(struct Cell ***,struct Sim *, struct MPIsetup *);
void cell_single_init_shear(struct Cell *, struct Sim *,int ,int ,int );
void cell_init_vortex(struct Cell ***,struct Sim *, struct MPIsetup *);
void cell_single_init_vortex(struct Cell *, struct Sim *,int ,int ,int );
void (*cell_init_ptr(struct Sim * ))(struct Cell *** , struct Sim * ,struct MPIsetup *);
void (*cell_single_init_ptr(struct Sim * ))(struct Cell * , struct Sim *,int,int,int );
///retrieve data
double cell_prim(struct Cell *, int);
double cell_grad(struct Cell *, int);
double cell_gradp(struct Cell *, int);
struct Cell *cell_single(struct Cell ***,int,int,int);
double cell_tiph(struct Cell *);
double cell_dphi(struct Cell *);
double cell_wiph(struct Cell *);
//modify cell data
void cell_add_cons(struct Cell *, int, double);
void cell_add_divB(struct Cell *, double);
void cell_add_GradPsi(struct Cell *, int, double);
void cell_add_wiph(struct Cell *, double);
void cell_add_src( struct Cell *** ,struct Sim * , struct GravMass * , double );
void cell_add_visc_src( struct Cell *** ,struct Sim * , double );
void cell_mult_psi(struct Cell *, double);
void cell_clean_pi(struct Cell *** ,struct Sim *);
void cell_update_phi( struct Cell *** , struct Sim * , double , double );
void cell_update_dphi( struct Cell *** ,struct Sim * );
//clear
void cell_clear_w(struct Cell ***,struct Sim * );
void cell_clear_divB( struct Cell ***,struct Sim * );
void cell_clear_GradPsi( struct Cell ***,struct Sim * );
//set w
//void cell_set_wcell(struct Cell ***,struct Sim *);
//void cell_set_wrigid(struct Cell ***,struct Sim *);
void cell_set(struct Cell ***,struct Sim *);
//processor syncs
void cell_syncproc_r( struct Cell *** , struct Sim *,struct MPIsetup *);
void cell_syncproc_z( struct Cell *** , struct Sim *,struct MPIsetup *);
//PLM
void cell_plm_rz( struct Cell *** ,struct Sim *, struct Face * , struct TimeStep * , int );
void cell_plm_p( struct Cell *** ,struct Sim * );
//boundary conditions
void cell_boundary_outflow_r( struct Cell *** , struct Face * ,struct Sim * ,struct MPIsetup *, struct TimeStep * );
void cell_boundary_outflow_z( struct Cell *** , struct Face * , struct Sim * ,struct MPIsetup *,struct TimeStep * );
void cell_boundary_fixed_r( struct Cell ***, struct Sim *,struct MPIsetup *,void (*)(struct Cell *,struct Sim *,int,int,int));
void cell_boundary_fixed_z( struct Cell ***, struct Sim *,struct MPIsetup *,void (*)(struct Cell *,struct Sim *,int,int,int));
void cell_bc_damp( struct Cell *** , struct Sim * , double ,void (*)(struct Cell *,struct Sim *,int,int,int));
//primitive-conservative conversion routines
void cell_calc_prim( struct Cell ***,struct Sim *);
void cell_prim2cons( double * , double * , double , double ,struct Sim *);
void cell_calc_cons(struct Cell ***,struct Sim *); 
void cell_cons2prim( double * , double * , double , double ,struct Sim *);
//miscellaneous
double cell_mindt( struct Cell *** , struct Sim * );
void cell_copy(struct Cell ***,struct Sim * );
void cell_adjust_RK_cons( struct Cell *** , struct Sim * , double );
void cell_set_prim(struct Cell ***,int,int,int,int,double);//this will morph into a checkpoint restart routine
void cell_set_tiph(struct Cell ***,int,int,int,double);//this will morph into a checkpoint restart routine
#endif
