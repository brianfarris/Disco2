#ifndef DIAGNOSTICS_H
#define DIAGNOSTICS_H
struct Diagnostics;
struct Grid;
struct Cell;
struct TimeStep;
struct MPIsetup;
#ifdef DIAGNOSTICS_PRIVATE_DEFS
struct Diagnostics{
  double **VectorDiag;
  double *ScalarDiag;
  double dtout;
  double toutprev;
  double toutprev_dump;
  int NUM_DIAG;
};
#endif
//create and destroy
struct Diagnostics *diagnostics_create(struct Grid *, struct TimeStep *);
void diagnostics_destroy(struct Diagnostics *,struct Grid *);
// set
void diagnostics_set(struct Diagnostics * theDiagnostics,struct Cell ***,struct Grid *,struct TimeStep *);
// print
void diagnostics_print(struct Diagnostics * ,struct TimeStep * ,struct Grid * ,struct MPIsetup * );
#endif
