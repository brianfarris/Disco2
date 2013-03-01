#ifndef RIEMANN_H
#define RIEMANN_H
struct Riemann;
struct Grid;
struct Cell;
struct Face;

#ifdef RIEMANN_PRIVATE_DEFS
struct Riemann {
  double *primL;
  double *primR;
  double *Uk;
  double *Ustar;
  double *F;
  double Sl;
  double Sr;
  double Ss;
  int state;
};
#endif

//create and destroy
struct Riemann *riemann_create(struct Grid * );
void riemann_destroy(struct Riemann *); 
//other routines
void riemann_set_primL(struct Riemann * theRiemann,int ,double );
void riemann_set_primR(struct Riemann * theRiemann,int ,double );
void riemann_set_vel(struct Riemann * ,double *,double ,double *,double ,double );
void riemann_set_state(struct Riemann * theRiemann,int );
void riemann_set_Ustar(struct Riemann * ,double *,double ,double *,double );
double * riemann_Uk(struct Riemann *);
double * riemann_Ustar(struct Riemann *);
void riemann_addto_flux_general(struct Riemann * ,double ,int );
void riemann_set_flux(struct Riemann * , double , double * ,double ,double );
double * riemann_prim(struct Riemann *,int );
int riemann_state(struct Riemann * );
double *riemann_F(struct Riemann * );
double riemann_Ss(struct Riemann * );
void riemann_blah( struct Cell * , struct Cell * ,struct Grid *, double , double , double ,double,double,double,double,int);
#endif
