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
  double Bk_face;
  double Psi_face;
};
#endif

//create and destroy
struct Riemann *riemann_create(struct Sim * );
void riemann_destroy(struct Riemann *); 
//other routines
void riemann_setup_rz(struct Riemann *,struct Face * , struct Sim *,int,int );
void riemann_setup_p(struct Riemann * ,struct Cell *** ,struct Sim * ,int ,int ,int,int );
void riemann_AddFlux(struct Riemann *, struct Sim *,double);
void riemann_add_to_cells(struct Riemann * ,struct Sim *,double );
#endif
