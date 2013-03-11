#ifndef TIMESTEP_H
#define TIMESTEP_H
struct TimeStep;
struct Cell;
struct Sim;
struct GravMass;
struct MPIsetup;

#ifdef TIMESTEP_PRIVATE_DEFS
struct TimeStep {
  double t;
  double dt;
  double RK;
  int *nri;
  int *nzk;
  int Nfr;
  int Nfz;
 };
#endif
// timestepping algorithms
void timestep_rk2(struct TimeStep * , struct Sim * , struct Cell *** , struct GravMass * , struct MPIsetup * );
//create and destroy
struct TimeStep *timestep_create(struct Sim * );
void timestep_destroy(struct TimeStep *); 
//adjust t, dt, and RK
void timestep_set_dt(struct TimeStep * , struct Cell *** , struct Sim * );
void timestep_set_t(struct TimeStep * ,double );
void timestep_update_t(struct TimeStep *);
void timestep_set_RK(struct TimeStep * ,double);
//take a substep
void timestep_substep(struct TimeStep * , struct Cell *** ,struct Sim * ,struct GravMass * ,struct MPIsetup *,double);
void timestep_update_Psi( struct TimeStep * , struct Cell *** , struct Sim *,struct MPIsetup * );
//access data 
double timestep_get_t(struct TimeStep *);
double timestep_dt(struct TimeStep * );
int timestep_nri(struct TimeStep *,int);
int timestep_nzk(struct TimeStep *,int);
//int timestep_Nfr(struct TimeStep *);
//int timestep_Nfz(struct TimeStep *);
//set data
//void timestep_set_Nfr(struct TimeStep *,struct Sim *);
//void timestep_set_Nfz(struct TimeStep *,struct Sim *);
void timestep_set_nri(struct TimeStep * ,int ,int );
void timestep_set_nzk(struct TimeStep * ,int ,int );
#endif
