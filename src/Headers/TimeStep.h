#ifndef TIMESTEP_H
#define TIMESTEP_H
struct TimeStep;
struct Cell;
struct Grid;
struct GravMass;
struct MPIsetup;

#ifdef TIMESTEP_PRIVATE_DEFS
struct TimeStep {
  double t;
  double dt;
  double RK;
  double T_MAX;
  int NUM_CHECKPOINTS;
};
#endif

//create and destroy
struct TimeStep *timestep_create();
void timestep_destroy(struct TimeStep *); 
//adjust t, dt, and RK
void timestep_set_dt(struct TimeStep * , struct Cell *** , struct Grid * );
void timestep_update_t(struct TimeStep *);
void timestep_set_RK(struct TimeStep * ,double);
//take a substep
void timestep_substep(struct TimeStep * , struct Cell *** ,struct Grid * ,struct GravMass * ,struct MPIsetup *,double);
void timestep_update_Psi( struct TimeStep * , struct Cell *** , struct Grid *,struct MPIsetup * );
//access data 
double timestep_get_t(struct TimeStep *);
double timestep_NUM_CHECKPOINTS(struct TimeStep *);
#endif
