#ifndef FLUX_H
#define FLUX_H
struct Flux;
struct MPIsetup;

#ifdef FLUX_PRIVATE_DEFS
struct Flux {
  double *primL;
  double *primR;
  double *Uk;
  double *Ustar;
  double *F;
  double Sl;
  double Sr;
  double Ss;
};
#endif

//create and destroy
struct Flux *flux_create(struct Grid * );
void flux_destroy(struct Flux *); 
//other routines
void flux_set_primL(struct Flux * theFlux,int q,double input);
void flux_set_primR(struct Flux * theFlux,int q,double input);
void flux_set_vel(struct Flux * ,double *,double ,double *,double ,double );
void flux_set_Ustar(struct Flux * ,double *,double ,double *,double ,double);
void flux_set_Uk(struct Flux *,double *);
double * flux_primL(struct Flux * );
double * flux_primR(struct Flux * );

