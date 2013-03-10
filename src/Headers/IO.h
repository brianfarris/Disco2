#ifndef IO_H
#define IO_H
struct IO;
struct Cell;
struct Sim;
#ifdef IO_PRIVATE_DEFS
struct IO{
  double **primitives;
};
#endif
//create and destroy
struct IO *io_create(struct Sim *);
void io_destroy(struct IO *);
//move data between theCells and IO buffer
void io_flattened_prim(struct IO *,struct Cell ***,struct Sim *);
void io_unflattened_prim(struct IO *,struct Cell ***,struct Sim *);
//calls to hdf5 routines
void io_hdf5_out(struct IO *,struct Sim *,char *);
void io_hdf5_in(struct IO *,struct Sim *,char * );
#endif
