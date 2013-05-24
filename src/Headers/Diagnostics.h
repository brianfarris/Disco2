#ifndef DIAGNOSTICS_H
#define DIAGNOSTICS_H
struct Diagnostics;
struct Sim;
struct Cell;
struct TimeStep;
struct MPIsetup;
#ifdef DIAGNOSTICS_PRIVATE_DEFS
struct Diagnostics{
  int * N_p_global;
  double **EquatDiag;
  double **VectorDiag;
  double *ScalarDiag;
  int offset_eq;
  int N_eq_cells;
  double dtout;
  double toutprev;
  double toutprev_dump;
  double tdiag_measure;
  double dtdiag_measure;
  double tdiag_dump;
  double dtdiag_dump;
  int NUM_DIAG;
};
#endif
//create and destroy
struct Diagnostics *diagnostics_create(struct Sim *, struct TimeStep *, struct MPIsetup *);
void diagnostics_destroy(struct Diagnostics *,struct Sim *);
// set
void diagnostics_reset(struct Diagnostics * theDiagnostics,struct Cell ***,struct Sim *,struct TimeStep *);
void diagnostics_set(struct Diagnostics * theDiagnostics,struct Cell ***,struct Sim *,struct TimeStep *,struct MPIsetup *,struct GravMass *);
// print
void diagnostics_print(struct Diagnostics * ,struct TimeStep * ,struct Sim * ,struct MPIsetup * );
// access
double diagnostics_tdiag_measure(struct Diagnostics * ); 
double diagnostics_tdiag_dump(struct Diagnostics * ); 
#endif
