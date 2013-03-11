#define DIAGNOSTICS_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include "../Headers/Diagnostics.h"
#include "../Headers/header.h"

double diagnostics_tdiag_measure(struct Diagnostics * theDiagnostics) {
  return(theDiagnostics->tdiag_measure);
}

double diagnostics_tdiag_dump(struct Diagnostics * theDiagnostics) {
  return(theDiagnostics->tdiag_dump);
}


