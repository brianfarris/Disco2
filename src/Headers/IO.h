#ifndef IO_H
#define IO_H
struct IO;
struct Cell;
struct TimeStep;
struct Sim;
struct GravMass;
#ifdef IO_PRIVATE_DEFS
struct IO{
  double **buffer;
  double tcheck;
  double dtcheck;
  int nfile;
  char filename[256];
  double GravMassBuffer[2][4];
};
#endif
//create and destroy
struct IO *io_create(struct Sim *);
void io_destroy(struct IO *);
//move data between theCells and IO buffer
void io_allocbuf(struct IO *,struct Sim *);
void io_deallocbuf(struct IO *);
void io_setbuf(struct IO *,struct Cell ***,struct Sim *,struct GravMass *);
void io_readbuf(struct IO *,struct Cell ***,struct Sim *,struct GravMass *);
//calls to hdf5 routines
void io_hdf5_out(struct IO *,struct Sim *,struct TimeStep *);
void io_hdf5_in(struct IO *,struct Sim *,struct TimeStep *);
// access
double io_tcheck(struct IO * ); 
int io_nfile(struct IO * ); 
//setup
void io_setup(struct IO * ,struct Sim * ,struct TimeStep * );
#endif
