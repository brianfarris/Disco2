#define DIAGNOSTICS_PRIVATE_DEFS
#include "../Headers/Diagnostics.h"
#include "../Headers/header.h"

struct Diagnostics *grid_create(struct Grid * theGrid) {
  struct Diagnostics * theDiagnostics = (struct Diagnostics *) malloc(sizeof(struct Diagnostics));
  theDiagnostics->VectorDiag =; 
  return theDiagnostics;
}

void grid_destroy(struct Diagnostics *theDiagnostics) {
  free(theDiagnostics->r_faces);
  free(theDiagnostics->z_faces);
  free(theDiagnostics->N_p);
  free(theDiagnostics);
}


