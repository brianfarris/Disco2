#define DIAGNOSTICS_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include "../Headers/Diagnostics.h"
#include "../Headers/TimeStep.h"
#include "../Headers/Sim.h"
#include "../Headers/MPIsetup.h"
#include "../Headers/header.h"

void diagnostics_print(struct Diagnostics * theDiagnostics,struct TimeStep * theTimeStep,struct Sim * theSim,struct MPIsetup * theMPIsetup){
  if (timestep_get_t(theTimeStep)>diagnostics_tdiag_dump(theDiagnostics)){
    double dt_dump = timestep_get_t(theTimeStep) - theDiagnostics->toutprev_dump;
    if(mpisetup_MyProc(theMPIsetup)==0){
      char DiagVectorFilename[256];
      char DiagScalarFilename[256];
      sprintf(DiagVectorFilename,"DiagVector_%6.4f.dat",timestep_get_t(theTimeStep));
      sprintf(DiagScalarFilename,"DiagScalar.dat");
      FILE * DiagVectorFile = fopen(DiagVectorFilename,"w");
      FILE * DiagScalarFile = fopen(DiagScalarFilename,"a");
      int i,n;
      for (i=0;i<sim_N_global(theSim,R_DIR);++i){
        for (n=0;n<theDiagnostics->NUM_DIAG+1;++n){
          fprintf(DiagVectorFile,"%e ",theDiagnostics->VectorDiag[i][n]/dt_dump);       
        }
        fprintf(DiagVectorFile,"\n");
      }
      fprintf(DiagScalarFile,"%e ",timestep_get_t(theTimeStep));       
      for (n=0;n<theDiagnostics->NUM_DIAG;++n){
        fprintf(DiagScalarFile,"%e ",theDiagnostics->ScalarDiag[n]/dt_dump);       
      }
      fprintf(DiagScalarFile,"\n");

      fclose(DiagVectorFile);
      fclose(DiagScalarFile);
    }
    //update output time;
    theDiagnostics->toutprev_dump = timestep_get_t(theTimeStep);

    //reset 
    int i,n;
    for (i=0;i<sim_N_global(theSim,R_DIR);++i){
      for (n=0;n<(1+theDiagnostics->NUM_DIAG);++n){
        theDiagnostics->VectorDiag[i][n]=0.0;
      }
    }
    for (n=0;n<theDiagnostics->NUM_DIAG;++n){
      theDiagnostics->ScalarDiag[n]=0.0;
    }

    theDiagnostics->tdiag_dump += theDiagnostics->dtdiag_dump;

  }
}
