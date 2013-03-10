#ifndef RIEMANN_H
#define RIEMANN_H
struct Riemann;
struct Sim;
struct Cell;
struct Face;

#ifdef RIEMANN_PRIVATE_DEFS
struct Riemann {
  struct Cell * cL;
  struct Cell * cR;
  double r;
  double dA;
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
struct Riemann *riemann_create(struct Sim * );
void riemann_destroy(struct Riemann *); 
//other routines
void riemann_set_vel(struct Riemann * ,struct Sim *,double *,double ,double *,double ,double );
void riemann_set_Ustar(struct Riemann * ,struct Sim *,double *,double ,double *,double );
void riemann_addto_flux_general(struct Riemann * ,double ,int );
void riemann_visc_flux(struct Riemann * ,struct Sim * , double * );
void riemann_set_flux(struct Riemann * ,struct Sim *, double , double * ,double ,double );
void riemann_setup_rz(struct Riemann *,struct Face * , struct Sim *,int );
void riemann_setup_p(struct Riemann * ,struct Cell *** ,struct Sim * ,int ,int ,int ,int );
//void riemann_hllc(struct Riemann *, struct Cell * , struct Cell * ,struct Sim *, double , double , double ,double,double,double,double,int);
void riemann_hllc(struct Riemann *, struct Sim *,double,int);
#endif
