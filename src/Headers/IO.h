#ifndef IO_H
#define IO_H
struct io;
struct cell;
struct grid;
#ifdef IO_PRIVATE_DEFS
struct io{
  double **primitives;
};
#endif
struct io *io_create(struct Grid *);
void io_destroy(struct io *,struct Grid *);
void io_flattened_prim(struct io *,struct Cell ***,struct Grid *);
void io_unflattened_prim(struct io *,struct Cell ***,struct Grid *);
void io_hdf5_out(struct io *,struct Grid *,char *);
void io_hdf5_in(struct io *,struct Grid * );
#endif
