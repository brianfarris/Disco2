#ifndef FACE_H
#define FACE_H
struct Face;
struct Grid;
struct Cell;

#ifdef FACE_PRIVATE_DEFS
struct Face{
  struct Cell * L;
  struct Cell * R;
  double r;
  double dphi;
  double dA;
  double deltaL;
  double deltaR;
  double cm;
};
#endif

struct Face *face_pointer( struct Face *,int);
struct Cell *face_L_pointer( struct Face *,int);
struct Cell *face_R_pointer( struct Face *,int);
double face_deltaL( struct Face *);
double face_deltaR( struct Face *);
double face_cm( struct Face *);
double face_dA( struct Face *);
struct Face *face_create_r(struct Cell *** ,struct Grid *,int *, int *);
struct Face *face_create_z(struct Cell *** ,struct Grid *,int *, int *);
void face_destroy(struct Face *);
void face_riemann_r( struct Face * , double );
void face_riemann_z( struct Face * , double );
#endif
