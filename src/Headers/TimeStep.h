#ifndef TIMESTEP_H
#define TIMESTEP_H
struct TimeStep;

#ifdef TIMESTEP_PRIVATE_DEFS
struct TimeStep {
  double t;
  double dt;
  double RK;
};
#endif

struct TimeStep *timestep_create();
void timestep_destroy(struct TimeStep *); 
#endif
