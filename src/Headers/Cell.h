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
  double * gradr;
  double * gradp;
  double * gradz;
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
void cell_init_vortex(struct Cell ***,struct Sim *, struct MPIsetup *);
void cell_single_init_vortex(struct Cell *, struct Sim *,int ,int ,int );
void cell_init_torus(struct Cell ***,struct Sim *, struct MPIsetup *);
void cell_single_init_torus(struct Cell *, struct Sim *,int ,int ,int );
void cell_init_bondi(struct Cell ***,struct Sim *, struct MPIsetup *);
void cell_single_init_bondi(struct Cell *, struct Sim *,int ,int ,int );
void cell_init_uniform(struct Cell ***,struct Sim *, struct MPIsetup *);
void cell_single_init_uniform(struct Cell *, struct Sim *,int ,int ,int );
void cell_init_shock1(struct Cell ***,struct Sim *, struct MPIsetup *);
void cell_single_init_shock1(struct Cell *, struct Sim *,int ,int ,int );
void cell_init_shock2(struct Cell ***,struct Sim *, struct MPIsetup *);
void cell_single_init_shock2(struct Cell *, struct Sim *,int ,int ,int );
void cell_init_shock3(struct Cell ***,struct Sim *, struct MPIsetup *);
void cell_single_init_shock3(struct Cell *, struct Sim *,int ,int ,int );
void cell_init_shock4(struct Cell ***,struct Sim *, struct MPIsetup *);
void cell_single_init_shock4(struct Cell *, struct Sim *,int ,int ,int );
void cell_init_isentrope(struct Cell ***,struct Sim *, struct MPIsetup *);
void cell_single_init_isentrope(struct Cell *, struct Sim *,int ,int ,int );
void cell_init_gbondi(struct Cell ***,struct Sim *, struct MPIsetup *);
void cell_single_init_gbondi(struct Cell *, struct Sim *,int ,int ,int );
void cell_init_gbondi2(struct Cell ***,struct Sim *, struct MPIsetup *);
void cell_single_init_gbondi2(struct Cell *, struct Sim *,int ,int ,int );
void cell_init_equil1(struct Cell ***,struct Sim *, struct MPIsetup *);
void cell_single_init_equil1(struct Cell *, struct Sim *,int ,int ,int );
void cell_init_equil2(struct Cell ***,struct Sim *, struct MPIsetup *);
void cell_single_init_equil2(struct Cell *, struct Sim *,int ,int ,int );
void cell_init_ssdisc(struct Cell ***,struct Sim *, struct MPIsetup *);
void cell_single_init_ssdisc(struct Cell *, struct Sim *,int ,int ,int );
void cell_init_ntdisc(struct Cell ***,struct Sim *, struct MPIsetup *);
void cell_single_init_ntdisc(struct Cell *, struct Sim *,int ,int ,int );
void cell_init_cnstdisc(struct Cell ***,struct Sim *, struct MPIsetup *);
void cell_single_init_cnstdisc(struct Cell *, struct Sim *,int ,int ,int );
void cell_init_cartshear(struct Cell ***,struct Sim *, struct MPIsetup *);
void cell_single_init_cartshear(struct Cell *, struct Sim *,int ,int ,int );
void cell_init_disctest(struct Cell ***,struct Sim *, struct MPIsetup *);
void cell_single_init_disctest(struct Cell *, struct Sim *,int ,int ,int );
void cell_init_accdisc(struct Cell ***,struct Sim *, struct MPIsetup *);
void cell_single_init_accdisc(struct Cell *, struct Sim *,int ,int ,int );
void cell_init_adaf(struct Cell ***,struct Sim *, struct MPIsetup *);
void cell_single_init_adaf(struct Cell *, struct Sim *,int ,int ,int );
void cell_init_adaf2(struct Cell ***,struct Sim *, struct MPIsetup *);
void cell_single_init_adaf2(struct Cell *, struct Sim *,int ,int ,int );
void cell_init_atmo(struct Cell ***,struct Sim *, struct MPIsetup *);
void cell_single_init_atmo(struct Cell *, struct Sim *,int ,int ,int );
void cell_init_orbit(struct Cell ***,struct Sim *, struct MPIsetup *);
void cell_single_init_orbit(struct Cell *, struct Sim *,int ,int ,int );
void cell_init_minidisc(struct Cell ***,struct Sim *, struct MPIsetup *);
void cell_single_init_minidisc(struct Cell *, struct Sim *,int ,int ,int );
void (*cell_init_ptr(struct Sim * ))(struct Cell *** , struct Sim * ,struct MPIsetup *);
void (*cell_single_init_ptr(struct Sim * ))(struct Cell * , struct Sim *,int,int,int );
///retrieve data
double cell_prim(struct Cell *, int);
double cell_gradr(struct Cell *, int);
double cell_gradz(struct Cell *, int);
double cell_gradp(struct Cell *, int);
struct Cell *cell_single(struct Cell ***,int,int,int);
double cell_tiph(struct Cell *);
double cell_dphi(struct Cell *);
double cell_wiph(struct Cell *);
//modify cell data
void cell_add_cons(struct Cell *, int, double);
void cell_add_wiph(struct Cell *, double);
void cell_clean_pi(struct Cell *** ,struct Sim *);
void cell_update_phi( struct Cell *** , struct Sim * , double , double );
void cell_update_dphi( struct Cell *** ,struct Sim * );
void (*cell_add_src)( struct Cell *** ,struct Sim * , struct GravMass * , double );
void cell_add_src_newt( struct Cell *** ,struct Sim * , struct GravMass * , double );
void cell_add_src_gr( struct Cell *** ,struct Sim * , struct GravMass * , double );
void cell_add_src_grdisc( struct Cell *** ,struct Sim * , struct GravMass * , double );
//clear
void cell_clear_w(struct Cell ***,struct Sim * );
//set w
void cell_set_w(struct Cell ***,struct Sim *);
//processor syncs
void cell_syncproc_r( struct Cell *** , struct Sim *,struct MPIsetup *);
void cell_syncproc_z( struct Cell *** , struct Sim *,struct MPIsetup *);
//PLM
void cell_plm_rz( struct Cell *** ,struct Sim *, struct Face * , struct TimeStep * , struct MPIsetup *, int );
void cell_plm_p( struct Cell *** ,struct Sim * );
//boundary conditions
void cell_boundary_outflow_r_inner( struct Cell *** , struct Face * ,struct Sim * ,struct MPIsetup *, struct TimeStep * );
void cell_boundary_outflow_r_outer( struct Cell *** , struct Face * ,struct Sim * ,struct MPIsetup *, struct TimeStep * );
void cell_boundary_outflow_z_bot( struct Cell *** , struct Face * , struct Sim * ,struct MPIsetup *,struct TimeStep * );
void cell_boundary_outflow_z_top( struct Cell *** , struct Face * , struct Sim * ,struct MPIsetup *,struct TimeStep * );
void cell_boundary_fixed_r_inner( struct Cell ***, struct Sim *,struct MPIsetup *,void (*)(struct Cell *,struct Sim *,int,int,int));
void cell_boundary_fixed_r_outer( struct Cell ***, struct Sim *,struct MPIsetup *,void (*)(struct Cell *,struct Sim *,int,int,int));
void cell_boundary_fixed_z_bot( struct Cell ***, struct Sim *,struct MPIsetup *,void (*)(struct Cell *,struct Sim *,int,int,int));
void cell_boundary_fixed_z_top( struct Cell ***, struct Sim *,struct MPIsetup *,void (*)(struct Cell *,struct Sim *,int,int,int));
void cell_boundary_ssprofile_r_inner( struct Cell ***, struct Sim *,struct MPIsetup *);
void cell_boundary_ssprofile_r_outer( struct Cell ***, struct Sim *,struct MPIsetup *);
void cell_boundary_linear_r_inner( struct Cell *** , struct Face * ,struct Sim * ,struct MPIsetup *, struct TimeStep * );
void cell_boundary_linear_r_outer( struct Cell *** , struct Face * ,struct Sim * ,struct MPIsetup *, struct TimeStep * );
void cell_boundary_nozzle( struct Cell ***, struct Sim *,struct MPIsetup *);
void cell_bc_damp( struct Cell *** , struct Sim * , double ,void (*)(struct Cell *,struct Sim *,int,int,int));
//primitive-conservative conversion routines
void cell_calc_prim( struct Cell ***,struct Sim *);
void cell_calc_cons(struct Cell ***,struct Sim *); 
void (*cell_prim2cons)( double * , double * , double * , double ,struct Sim *);
void (*cell_cons2prim)( double * , double * , double * , double ,struct Sim *);
void cell_prim2cons_newt( double * , double * , double * , double ,struct Sim *);
void cell_cons2prim_newt( double * , double * , double * , double ,struct Sim *);
void cell_prim2cons_gr( double * , double * , double * , double ,struct Sim *);
void cell_cons2prim_gr( double * , double * , double * , double ,struct Sim *);
void cell_prim2cons_grdisc( double * , double * , double * , double ,struct Sim *);
void cell_cons2prim_grdisc( double * , double * , double * , double ,struct Sim *);
//miscellaneous
double (*cell_mindt)( struct Cell *** , struct Sim * );
double cell_mindt_newt( struct Cell *** , struct Sim * );
double cell_mindt_gr( struct Cell *** , struct Sim * );
double cell_mindt_grvisc1( struct Cell *** , struct Sim * );
double cell_mindt_grdisc( struct Cell *** , struct Sim * );
void cell_copy(struct Cell ***,struct Sim * );
void cell_adjust_RK_cons( struct Cell *** , struct Sim * , double );
void cell_set_prim(struct Cell ***,int,int,int,int,double);//this will morph into a checkpoint restart routine
void cell_set_tiph(struct Cell ***,int,int,int,double);//this will morph into a checkpoint restart routine
void cell_print(struct Cell *** ,int ,int ,int );
void cell_print_all(struct Cell *** , struct Sim *);
void cell_print_nan(struct Cell ***theCells, char label[], struct Sim *theSim);
void cell_metric_horizon_jazz(struct Cell ***, struct Sim *);
#endif
