#ifndef IO_H
#define IO_H
struct IO;
struct cell;
struct grid;
#ifdef IO_PRIVATE_DEFS
struct IO{
  double **primitives;
};
#endif
//create and destroy
struct IO *io_create(struct Grid *);
void io_destroy(struct IO *,struct Grid *);
//move data between theCells and IO buffer
void io_flattened_prim(struct IO *,struct Cell ***,struct Grid *);
void io_unflattened_prim(struct IO *,struct Cell ***,struct Grid *);
//calls to hdf5 routines
void io_hdf5_out(struct IO *,struct Grid *,char *);
void io_hdf5_in(struct IO *,struct Grid *,char * );
#endif
