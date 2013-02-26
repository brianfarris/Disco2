#ifndef GRID_H
#define GRID_H
struct Grid;
struct MPIsetup;

#ifdef GRID_PRIVATE_DEFS
struct Grid {
  double *r_faces;
  double *z_faces;
  int *N_p;
  int N_r;
  int N_z;
  int Nghost_rmin;
  int Nghost_rmax;
  int Nghost_zmin;
  int Nghost_zmax;
  int Ncells;
  int Ncells_global;
  int offset;
  int N_r_global;
  int N_z_global;
  int ng;
  double RMIN;
  double RMAX;
  double ZMIN;
  double ZMAX;

};
#endif

//create and destroy
struct Grid *grid_create(struct MPIsetup * );
void grid_destroy(struct Grid *); 
//access grid data
int grid_N_p(struct Grid *,int);
double grid_r_faces(struct Grid *,int);
double grid_z_faces(struct Grid *,int);
int grid_N_r(struct Grid *);
int grid_N_z(struct Grid *);
int grid_Ncells(struct Grid *);
int grid_Ncells_global(struct Grid *);
int grid_offset(struct Grid *);
int grid_Nghost_rmin(struct Grid *);
int grid_Nghost_rmax(struct Grid *);
int grid_Nghost_zmin(struct Grid *);
int grid_Nghost_zmax(struct Grid *);
//set grid data
void grid_set_N_p(struct Grid *);
void grid_set_rz(struct Grid *,struct MPIsetup *);
void grid_set_Ncells_and_offset(struct Grid *,struct MPIsetup *);
#endif
