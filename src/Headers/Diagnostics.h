#ifndef DIAGNOSTICS_H
#define DIAGNOSTICS_H
#define NUM_DIAG 4
struct Diagnostics;
#ifdef DIAGNOSTICS_PRIVATE_DEFS
struct Diagnostics{
  double ScalarDiag[NUM_DIAG]; 
  double *VectorDiag[NUM_DIAG];
  double **PlanarDiag[NUM_DIAG];
};
#endif
//create and destroy
struct Diagnostics *diagnostics_create(struct Grid *);
void diagnostics_destroy(struct Diagnostics *);
#endif
