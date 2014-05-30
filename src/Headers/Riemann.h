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
  double pos[3];
  double r_cell_L;
  double r_cell_R;
  double x_cell_L;
  double x_cell_R;
  double dA;
  double cm;
  double *primL;
  double *primR;
  double *UL;
  double *UR;
  double *Ustar;
  double *FL;
  double *FR;
  double *Fstar;
  double *F;
  double Sl;
  double Sr;
  double Ss;
  int state;
  int n[3];
};
#endif

//create and destroy
struct Riemann *riemann_create(struct Sim * );
void riemann_destroy(struct Riemann *); 
//other routines
void riemann_setup_rz(struct Riemann *,struct Face * , struct Sim *,int,int );
void riemann_setup_p(struct Riemann * ,struct Cell *** ,struct Sim * ,int ,int ,int,int );
void riemann_AddFlux(struct Riemann *, struct Sim *,double);

void (*riemann_set_flux)(struct Riemann *, struct Sim *, double, int);
void riemann_set_flux_newt(struct Riemann *, struct Sim *, double, int);
void riemann_set_flux_gr(struct Riemann *, struct Sim *, double, int);

void (*riemann_set_vel)(struct Riemann *, struct Sim *, double, double);
void riemann_set_vel_newt(struct Riemann *, struct Sim *, double, double);
void riemann_set_vel_gr(struct Riemann *, struct Sim *, double, double);
  
void riemann_visc_flux(struct Riemann *, struct Sim *);
void riemann_visc_flux_LR(struct Riemann *, struct Sim *, int, double *);
#endif
