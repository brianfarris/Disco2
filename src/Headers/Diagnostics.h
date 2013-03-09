#ifndef DIAGNOSTICS_H
#define DIAGNOSTICS_H
struct Diagnostics;
struct Grid;
struct TimeStep;
#ifdef DIAGNOSTICS_PRIVATE_DEFS
struct Diagnostics{
  double **VectorDiag;
  double *ScalarDiag;
  double dtout;
  double toutprev;
};
#endif
//create and destroy
struct Diagnostics *diagnostics_create(struct Grid *, struct TimeStep *);
void diagnostics_destroy(struct Diagnostics *,struct Grid *);
// set
void diagnostics_set(struct Diagnostics * theDiagnostics,struct Cell ***,struct Grid *,struct TimeStep *);
#endif
