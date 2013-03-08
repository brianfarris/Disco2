#ifndef DIAGNOSTICS_H
#define DIAGNOSTICS_H
struct Diagnostics;
#ifdef DIAGNOSTICS_PRIVATE_DEFS
struct Diagnostics{
  double **VectorDiag;
  double *ScalarDiag; 
};
#endif
//create and destroy
struct Diagnostics *diagnostics_create(struct Grid *);
void diagnostics_destroy(struct Diagnostics *);
// set
void diagnostics_set(struct Diagnostics * theDiagnostics,struct Cell ***,struct Grid *);
#endif
