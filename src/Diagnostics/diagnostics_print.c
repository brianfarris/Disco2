#define DIAGNOSTICS_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include "../Headers/Diagnostics.h"
#include "../Headers/TimeStep.h"
#include "../Headers/Sim.h"
#include "../Headers/MPIsetup.h"
#include "../Headers/header.h"
#include "hdf5.h"


void diagnostics_print(struct Diagnostics * theDiagnostics,struct TimeStep * theTimeStep,struct Sim * theSim,struct MPIsetup * theMPIsetup){
  if (timestep_get_t(theTimeStep)>diagnostics_tdiag_dump(theDiagnostics)){
    double dt_dump = timestep_get_t(theTimeStep) - theDiagnostics->toutprev_dump;
    if(mpisetup_MyProc(theMPIsetup)==0){
      int NUM_EQ = theDiagnostics->NUM_DIAG+2;
      char DiagEquatFilename[256];
      char DiagVectorFilename[256];
      char DiagScalarFilename[256];
      sprintf(DiagEquatFilename,"DiagEquat_%08.4f.h5",timestep_get_t(theTimeStep));
      sprintf(DiagVectorFilename,"DiagVector_%08.4f.dat",timestep_get_t(theTimeStep));
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

      // HDF5 APIs definitions
      hid_t       file_id, dset_id;         /* file and dataset identifiers */
      hid_t       filespace, memspace;      /* file and memory dataspace identifiers */
      hsize_t     dimsf1[1];                 /* 1d dataset dimensions */
      hsize_t     dimsf2[2];                 /* 2d dataset dimensions */
      hid_t	plist_id;                 /* property list identifier */
      herr_t	status;
      // MPI variables
      MPI_Info info  = MPI_INFO_NULL;

      // Create a new file collectively and release property list identifier.
      file_id = H5Fcreate(DiagEquatFilename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT /*plist_id*/);

      // ***************************
      // Now we save the equatorial data
      // ***************************

      // Create the dataspace for the dataset.
      dimsf2[0] = theDiagnostics->N_eq_cells;
      dimsf2[1] = NUM_EQ;
      filespace = H5Screate_simple(2, dimsf2, NULL); 
      memspace  = H5Screate_simple(2, dimsf2, NULL); 


      dset_id = H5Dcreate1(file_id, "EQUAT", H5T_NATIVE_DOUBLE, filespace,H5P_DEFAULT);//,plist_id);
      H5Sclose(filespace);

      filespace = H5Dget_space(dset_id);

      status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT , &(theDiagnostics->EquatDiag[0][0]));

      // Close/release resources.
      H5Dclose(dset_id);
      H5Sclose(filespace);
      H5Sclose(memspace);
      H5Fclose(file_id);
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
