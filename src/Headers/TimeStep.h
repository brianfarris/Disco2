#ifndef TIMESTEP_H
#define TIMESTEP_H
struct TimeStep;
struct Cell;
struct Grid;
struct GravMass;

#ifdef TIMESTEP_PRIVATE_DEFS
struct TimeStep {
  double t;
  double dt;
  double RK;
};
#endif

struct TimeStep *timestep_create();
void timestep_destroy(struct TimeStep *); 
void timestep_set_dt(struct TimeStep * , struct Cell *** , struct Grid * );
void timestep_update_t(struct TimeStep *);
void timestep_set_RK(struct TimeStep * ,double);
void timestep_substep(struct TimeStep * , struct Cell *** ,struct Grid * ,struct GravMass * ,double);
void timestep_update_Psi( struct TimeStep * , struct Cell *** , struct Grid * );
double timestep_get_t(struct TimeStep *);
#endif
