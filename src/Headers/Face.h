#ifndef FACE_H
#define FACE_H
struct Face;
struct Sim;
struct Cell;
struct TimeStep;

#ifdef FACE_PRIVATE_DEFS
struct Face{
  struct Cell * L;
  struct Cell * R;
  double pos[3];
  double dphi;
  double dA;
  double deltaL;
  double deltaR;
  double cm;
};
#endif

//create and destroy
struct Face *face_create(struct Cell *** ,struct Sim *,struct TimeStep *, int);
void face_destroy(struct Face *);
//access face data
struct Cell *face_L_pointer( struct Face *,int);
struct Cell *face_R_pointer( struct Face *,int);
double face_deltaL( struct Face *,int);
double face_deltaR( struct Face *,int);
double face_cm( struct Face *,int);
double face_dA( struct Face *,int);
double face_r( struct Face *,int);
double face_phi( struct Face *,int);
double face_z( struct Face *,int);
#endif
